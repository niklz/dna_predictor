library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(stringr)
library(forcats)
library(ranger)
library(ggplot2)
library(vip)
library(tidymodels)
library(probably)
library(embed) # For step_lmer (High-cardinality encoding)
library(themis)
library(tictoc)
library(future)
library(doFuture)
library(GGally)
library(config)
library(openxlsx)
library(DBI)
library(odbc)
library(progressr)
library(readxl)

conf <- config::get()

# -------------------------------------------------------------------------
# CUSTOM HELPER FUNCTIONS
# -------------------------------------------------------------------------
collect_ethnicity <- function(x) {
  dplyr::case_when(
    x %in% c("Not Stated", "Not Known", "Not Collected At This Time", "Not Set") ~ "Unknown",
    x %in% c("White - British", "White - Any other White background", "White - Irish", "Any other White background") ~ "White",
    x %in% c("Asian or Asian British - Indian", "Asian - Indian or British Indian", 
             "Asian or Asian British - Pakistani", "Asian - Pakistani or British Pakistani",
             "Asian or Asian British - Bangladeshi", "Asian - Bangladeshi or British Bangladeshi",
             "Asian or Asian British - Any other Asian background", "Asian - Any other Asian background") ~ "Asian",
    x %in% c("Black or Black British - Caribbean", "Black - Caribbean or Black British Caribbean",
             "Black Caribbean", "Black or Black British - African", 
             "Black  - African or British African", "Black or Black British - Any other Black background",
             "Black - Any other Black background") ~ "Black",
    x %in% c("Mixed - White and Black African", "Mixed - Any other background", 
             "Mixed - White and Black Caribbean", "Mixed - White and Asian", "Any other Mixed Background") ~ "Mixed",
    x %in% c("Other Ethnic Group - Chinese", "Chinese", "Any Other Ethnic Group", "Any other Ethnic Group") ~ "Other",
    TRUE ~ "Unknown"
  )
}

aggregate_ethnicity_high_level <- function(x) {
  x_clean <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_clean %in% c("unknown", "not stated", "not known", "not collected at this time", "not set") ~ "Unknown",
    grepl("^white", x_clean) ~ "White",
    grepl("^asian", x_clean) ~ "Asian",
    grepl("^black", x_clean) ~ "Black",
    grepl("^mixed", x_clean) ~ "Mixed",
    grepl("chinese", x_clean) | grepl("other", x_clean) ~ "Chinese",
    TRUE ~ "Unknown"
  )
}

parse_appt_date <- function(x) {
  if (inherits(x, "POSIXt") || inherits(x, "Date")) {
    return(as.Date(x))
  }
  parsed <- lubridate::parse_date_time(x, orders = c("dmy HMS", "ymd HMS", "dmy", "ymd"))
  return(as.Date(parsed))
}

parse_appt_time <- function(time_val) {
  if (inherits(time_val, "POSIXt")) {
    return(format(time_val, "%H:%M"))
  }
  
  sapply(as.character(time_val), function(val) {
    if (is.na(val) || val == "" || val == "NA") return("")
    if (grepl(" ", val)) {
      parts <- strsplit(val, " ")[[1]]
      if (length(parts) >= 2 && grepl("^\\d{1,2}:\\d{2}", parts[2])) {
        return(substring(parts[2], 1, 5))
      }
    }
    if (grepl("^\\d{1,2}:\\d{2}", val)) {
      return(substring(val, 1, 5))
    }
    return("")
  }, USE.NAMES = FALSE)
}


apply_custom_feature_engineering <- function(data) {
  data %>%
    mutate(
      across(c(age_at_appointment, distance_km, lead_time_days), as.numeric),
      # --- 1. Your Original Feature Engineering ---
      ethnicity_clean = collect_ethnicity(ethnicity),
      ethnicity_group = aggregate_ethnicity_high_level(ethnicity_clean),
      # Set the baseline reference level to the most frequent category
      ethnicity_group = relevel(factor(ethnicity_group), ref = "White"),
      appt_date = parse_appt_date(appt_month),
      # appt_dow = factor(weekdays(appt_date)),
      appt_month_num = as.factor(format(appt_date, "%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      lead_time_days_log = log1p(pmax(0, lead_time_days)),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      appt_hour_sin = sin(2 * pi * appt_hour / 24),
      appt_hour_cos = cos(2 * pi * appt_hour / 24),
      has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
    ) %>%
    rename(imd = index_multiple_deprivation_decile) %>%
    mutate(
      
      # Coerce Character variables 
      across(c(local_spec_code, national_spec_code, imd), as.character),
      
      # Coerce Integer variables (protects ID and indicators)
      dim_patient_id = as.integer(dim_patient_id),
      nfa_ind        = as.integer(nfa_ind),
      appt_wknd_ind  = as.integer(appt_wknd_ind),
      
      # Bulk coerce all 20 clinical accessibility flags ('a_') to integers
      across(starts_with("a_"), as.integer)
    ) %>%
    # --- 3. Clean up raw columns (Replaces step_rm() to avoid workflow leaks) ---
    select(
      -appt_hour, -lead_time_days, -appt_date, -appt_month, 
      -prev_dna_ly, -ethnicity, -ethnicity_clean, -age_at_appointment
    )
}


# =====================================================================
# find_optimal_ml_thresholds()
# =====================================================================
#' Find operationally- and statistically-feasible ML risk thresholds
#'
#' Given a precision-recall curve and the historical distribution of model
#' scores, this function identifies two "safe zone" thresholds that bound
#' the proportion of the clinic population that should be flagged by the
#' model:
#'
#'   - a FLOOR threshold: the minimum proportion needed to reach your
#'     target trial sample size (n_req) within the trial window (w weeks)
#'   - a CEILING threshold: the maximum proportion the care team can
#'     operationally handle, given call capacity (c_max) per week
#'
#' Any threshold you choose between these two values should simultaneously
#' satisfy statistical power requirements and stay within staff capacity.
#'
#' @param pr_curve_data A tibble/data frame from yardstick::pr_curve(),
#'   containing (at minimum) the columns `.threshold` and `precision`.
#'   Typically produced via:
#'   collect_predictions(...) %>% pr_curve(truth = ..., .pred_DNA)
#'
#' @param historical_scores Numeric vector of historical ML risk scores
#'   (predicted probabilities) from your cohort. Used to empirically
#'   estimate, for each candidate threshold, what proportion of real
#'   patients would be flagged. Must be a plain numeric vector — if
#'   extracting from collect_predictions(), remember to pull() the score
#'   column first, e.g. collect_predictions(...) %>% pull(.pred_DNA).
#'
#' @param v_week Numeric. Expected weekly clinic volume (patients/week)
#'   eligible for enrollment. Default 308.
#'
#' @param w Numeric. Planned trial duration in weeks. Default 26.
#'
#' @param n_req Numeric. Target total sample size required for adequate
#'   statistical power. Default 800.
#'
#' @param a_int Numeric in (0, 1]. Allocation ratio to the intervention
#'   arm (e.g. 0.50 = 50% of enrolled patients receive the intervention).
#'   Default 0.50.
#'
#' @param r_call Numeric in (0, 1]. Expected non-response rate among
#'   intervention-arm patients that will require a follow-up phone call
#'   from a coordinator. Default 0.70.
#'
#' @param c_max Numeric. Maximum number of calls coordinators can
#'   realistically handle per week (operational capacity constraint).
#'   Default 30.
#'
#' @return An (invisible) list with:
#'   \describe{
#'     \item{floor_threshold}{Model score threshold (>=) that most closely
#'       achieves the statistical floor proportion.}
#'     \item{ceiling_threshold}{Model score threshold (>=) that most
#'       closely achieves the operational ceiling proportion.}
#'     \item{floor_prop}{Empirical proportion of historical_scores flagged
#'       at floor_threshold.}
#'     \item{ceiling_prop}{Empirical proportion of historical_scores
#'       flagged at ceiling_threshold.}
#'     \item{floor_target_prop}{The theoretical target proportion derived
#'       from n_req, v_week, and w (before matching to an actual
#'       threshold).}
#'     \item{ceiling_target_prop}{The theoretical target proportion
#'       derived from c_max, v_week, a_int, and r_call (before matching).}
#'     \item{feasible}{Logical. TRUE if floor_prop <= ceiling_prop, i.e.
#'       a valid "safe zone" of thresholds exists between the two bounds.}
#'     \item{enriched_pr_curve}{The full pr_curve_data with an added
#'       prop_flagged column, useful for downstream plotting or
#'       sensitivity analysis.}
#'   }
#'
#' @examples
#' \dontrun{
#' optimal_thresholds <- find_optimal_ml_thresholds(
#'   pr_curve_data     = collect_predictions(model_tune$tune_res,
#'                          parameters = model_tune$best_params) %>%
#'                        pr_curve(truth = dna_outcome, .pred_DNA),
#'   historical_scores = collect_predictions(model_tune$tune_res,
#'                          parameters = model_tune$best_params) %>%
#'                        pull(.pred_DNA),
#'   v_week = 308,
#'   w      = 26,
#'   n_req  = 800
#' )
#' }
find_optimal_ml_thresholds <- function(
    pr_curve_data,
    historical_scores,
    v_week = 308,
    w = 26,
    n_req = 800,
    a_int = 0.50,
    r_call = 0.70,
    c_max = 30
) {
  
  # ---------------------------------------------------------
  # 1. Calculate Operational & Statistical Bounds (Proportions)
  # ---------------------------------------------------------
  floor_prop   <- n_req / (v_week * w)
  ceiling_prop <- c_max / (v_week * a_int * r_call)
  is_feasible  <- floor_prop <= ceiling_prop
  
  cat("====================================================\n")
  cat("      ML THRESHOLD SAFE ZONE EXTRACTOR              \n")
  cat("====================================================\n")
  cat(sprintf("Statistical Floor Target:   %.2f%% of clinic (Min N = %d over %d weeks)\n",
              floor_prop * 100, n_req, w))
  cat(sprintf("Operational Ceiling Target: %.2f%% of clinic (Max Calls = %d/week)\n",
              ceiling_prop * 100, c_max))
  cat("----------------------------------------------------\n")
  
  if (!is_feasible) {
    cat("[!] WARNING: Statistical floor exceeds operational ceiling!\n")
    cat("    No single threshold can satisfy both constraints -\n")
    cat("    consider raising c_max, extending w, or lowering n_req.\n\n")
  }
  
  # ---------------------------------------------------------
  # 2. Enrich PR Curve with Empirical Proportion Flagged (VECTORIZED)
  # ---------------------------------------------------------
  sorted_scores <- sort(historical_scores)
  n_hist <- length(sorted_scores)
  
  enriched_pr <- pr_curve_data %>%
    filter(.threshold != Inf & !is.na(.threshold)) %>%
    mutate(
      n_less       = findInterval(.threshold, sorted_scores, left.open = TRUE),
      prop_flagged = (n_hist - n_less) / n_hist
    ) %>%
    select(-n_less)
  
  # ---------------------------------------------------------
  # 3. Match Target Proportions to Closest PR Curve Thresholds
  # ---------------------------------------------------------
  floor_match <- enriched_pr %>%
    slice_min(abs(prop_flagged - floor_prop), n = 1, with_ties = FALSE)
  
  ceiling_match <- enriched_pr %>%
    slice_min(abs(prop_flagged - ceiling_prop), n = 1, with_ties = FALSE)
  
  # Estimated weekly/total volumes at each matched threshold, for context
  floor_n_per_week   <- floor_match$prop_flagged * v_week
  ceiling_n_per_week <- ceiling_match$prop_flagged * v_week
  floor_total_n      <- floor_n_per_week * w
  ceiling_est_calls  <- ceiling_n_per_week * a_int * r_call
  
  # ---------------------------------------------------------
  # 4. Console Report
  # ---------------------------------------------------------
  cat("RECOMMENDED MODEL THRESHOLD RANGE:\n\n")
  cat(sprintf("-> Power Floor Threshold:      >= %.4f\n", floor_match$.threshold))
  cat(sprintf("   Flags top %.1f%% of clinic (~%.0f pts/week, ~%.0f over %d weeks) | Precision: %.3f\n",
              floor_match$prop_flagged * 100, floor_n_per_week, floor_total_n, w, floor_match$precision))
  cat(sprintf("-> Capacity Ceiling Threshold:  >= %.4f\n", ceiling_match$.threshold))
  cat(sprintf("   Flags top %.1f%% of clinic (~%.0f pts/week, ~%.0f est. calls/week) | Precision: %.3f\n",
              ceiling_match$prop_flagged * 100, ceiling_n_per_week, ceiling_est_calls, ceiling_match$precision))
  cat("----------------------------------------------------\n")
  cat(sprintf("Safe zone: %s\n", if (is_feasible) "VALID (floor <= ceiling)" else "INVALID (floor > ceiling)"))
  cat("====================================================\n\n")
  
  # ---------------------------------------------------------
  # 5. Return structured list for downstream pipeline scripts
  # ---------------------------------------------------------
  return(invisible(list(
    floor_threshold     = floor_match$.threshold,
    ceiling_threshold   = ceiling_match$.threshold,
    floor_prop          = floor_match$prop_flagged,
    ceiling_prop        = ceiling_match$prop_flagged,
    floor_target_prop   = floor_prop,
    ceiling_target_prop = ceiling_prop,
    feasible            = is_feasible,
    enriched_pr_curve   = enriched_pr
  )))
}


# Centralised recipe builder to prevent train-serve skew and environment leakage [cite: 11]
build_trial_recipe <- function(data_template, target_col, fct_other_prp = 0.05) {
  
  formula_obj <- as.formula(paste(target_col, "~ ."))
  # Crucial: Strip parent environment references to keep the object lightweight [cite: 11]
  environment(formula_obj) <- baseenv() 
  
  recipe(formula_obj, data = data_template) %>%
    update_role(dim_patient_id, new_role = "id") %>%
    step_novel(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors(), -imd) %>%
    step_other(all_nominal_predictors(), threshold = fct_other_prp) %>%
    
    # Target encode high-cardinality nominative features [cite: 236]
    step_lencode_mixed(
      c(
        "clinic_location",
        "clinic_code",
        "site_code",
        "registered_gp_practice",
        "national_spec_code"
      ),
      outcome = vars(!!sym(target_col))
    ) %>%
    
    step_nzv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors())
}

standardise_phone_number <- function(phone_col) {
  phone_char <- stringr::str_trim(as.character(phone_col))
  
  needs_zero <- !is.na(phone_char) & 
    phone_char != "" & 
    phone_char != "NA" & 
    !stringr::str_starts(phone_char, "0")
  
  output_phone <- phone_char
  output_phone[needs_zero] <- stringr::str_c("0", phone_char[needs_zero])
  output_phone[is.na(phone_char) | phone_char == "" | phone_char == "NA"] <- ""
  
  return(output_phone)
}


generate_appointment_manifest <- function(new_appointments, op_threshold, rf_model, rf_calibrator, conf) {
  
  # A. Force lowercase column names to prevent case mismatches [cite: 68]
  new_appointments <- rename_with(new_appointments, .fn = str_to_lower)
  
  # B. Standardise tumor site spelling variations dynamically [cite: 109, 121]
  if ("tumorsite" %in% colnames(new_appointments)) {
    new_appointments <- rename(new_appointments, tumoursite = tumorsite)
  }
  if (!"tumoursite" %in% colnames(new_appointments)) {
    new_appointments$tumoursite <- "Unknown" # Fallback group to prevent crashes
  }
  
  # B. Feature engineering
  message("Engineering features for incoming cohort...")
  engineered_appointments <- apply_custom_feature_engineering(new_appointments)
  
  # C. Prediction and calibration
  message("Generating calibrated risk predictions...")
  raw_predictions <- predict(rf_model, new_data = engineered_appointments, type = "prob")
  calibrated_predictions <- raw_predictions %>% cal_apply(rf_calibrator)
  
  # D. Stratify risk and allocate trial arms (Within-Clinic Block Randomisation)
  message("Stratifying patient risk profiles and allocating cohorts within clinical groups...")
  allocation_ratio <- ifelse(!is.null(conf$trial_allocation_ratio), 
                             conf$trial_allocation_ratio, 
                             0.50)
  
  final_manifest <- new_appointments %>%
    mutate(
      dna_probability = calibrated_predictions$.pred_DNA,
      risk_profile = ifelse(dna_probability >= op_threshold, "high-risk", "low-risk"),
      trial_arm = "not in trial"
    )
  
  # Group-based stratified randomisation at the tumoursite level
  final_manifest <- final_manifest %>%
    group_by(tumoursite) %>%
    group_modify(~ {
      high_risk_indices <- which(.x$risk_profile == "high-risk")
      n_high_risk       <- length(high_risk_indices)
      
      if (n_high_risk > 0) {
        n_control <- round(n_high_risk * allocation_ratio)
        n_interv  <- n_high_risk - n_control
        
        # Perfect 50/50 shuffle block
        balanced_arms <- sample(c(
          rep("control", n_control),
          rep("intervention", n_interv)
        ))
        
        .x$trial_arm[high_risk_indices] <- balanced_arms
      }
      .x
    }) %>%
    ungroup()
  
  # E. Generate primary keys, collapse flags, and map variables (31 columns total)
  message("Generating unique keys and collapsing vulnerability flags...")
  
  vulnerability_cols <- c(
    "a_ld", "a_autism", "a_interpreter_req_bsl", "a_interpreter_req_lang", "a_balance",  
    "a_cognitive_impairment", "a_mobility_restriction", "a_hear_vis_impaired", "a_dementia",  
    "a_depression", "a_downs_syndrome", "a_long_standing_condition", "a_makaton",  
    "a_mild_cognitive_impairment", "a_memory_impairment", "a_mood_disorder", "a_other_disability",  
    "a_psychosis", "a_severe_anxiety", "a_wheelchair_user"
  )
  
  vulnerability_labels <- c(
    "Learning disability", "Autism", "BSL interpreter required", "Language interpreter required",
    "Balance impairment", "Cognitive impairment", "Mobility restriction", "Sensory impairment",
    "Dementia", "Depression", "Down's syndrome", "Long-standing condition", "Makaton communication",
    "Mild cognitive impairment", "Memory impairment", "Mood disorder", "Other disability",
    "Psychosis", "Severe anxiety", "Wheelchair user"
  )
  
  final_manifest <- final_manifest %>%
    rowwise() %>%
    mutate(
      v_vals = list(c_across(all_of(vulnerability_cols))),
      accessibility_flags = paste(vulnerability_labels[which(v_vals == 1)], collapse = ", ")
    ) %>%
    ungroup() %>%
    select(-all_of(vulnerability_cols), -v_vals) %>%
    mutate(
      appointment_id          = op_appt_id,
      patient_id              = pasid, 
      nhs_number              = nhsnumber,
      
      # Clean clear Date column
      appointment_date        = if (inherits(appt_dttm, "POSIXt")) {
        format(appt_dttm, "%Y-%m-%d")
      } else {
        parsed_dt <- lubridate::parse_date_time(
          appt_dttm, 
          orders = c("dmy HM", "dmy HMS", "dmy")
        )
        format(parsed_dt, "%Y-%m-%d")
      },
      
      # New clear Time column
      appointment_time        = {
        # Detect if any columns contain date/time to find the second column
        time_cols <- grep("time", colnames(new_appointments), value = TRUE)
        time_cols <- grep("appt|date", time_cols, value = TRUE)
        time_cols <- time_cols[!time_cols %in% c("lead_time_days", "lead_time_days_log")]
        
        raw_vals <- if (length(time_cols) > 0) {
          new_appointments[[time_cols[1]]]
        } else {
          appt_dttm
        }
        parse_appt_time(raw_vals)
      },
      
      appointment_day_of_week = format(as.Date(appt_dttm), "%A"),
      clinic_code             = clinic_code,  
      tumoursite              = tumoursite,
      gp_practice             = registered_gp_practice,
      age                     = age_at_appointment,
      sex                     = gender,
      ethnicity               = ethnicity,
      imd                     = index_multiple_deprivation_decile,
      dna_risk                = dna_probability,
      risk_label              = risk_profile,
      
      date_model_run          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      model_version           = ifelse(!is.null(conf$model_ver), conf$model_ver, "v1.1"),
      
      patient_landline_number = standardise_phone_number(homephonenumber),
      patient_mobile_number   = standardise_phone_number(mobilephonenumber),
      
      tier1_text_timestamp    = "",
      tier2_text_timestamp    = "",
      t2_patient_intent       = "",
      tier3_phone_attempt     = "",
      tier3_call_timestamp    = "",
      number_of_call_attempts = "",
      call_duration_minutes   = "",
      t3_patient_intent       = "",
      call_notes              = "",
      appointment_outcome     = ""
    )
  
  return(final_manifest)
}


generate_excel_manifest <- function(manifest, org_name, op_threshold, output_dir = "outputs") {
  
  # A. Dynamically construct file paths
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  manifest_csv_path   <- file.path(output_dir, sprintf("%s_weekly_scored_appointments_%s.csv", org_name, current_date))
  action_list_path    <- file.path(output_dir, sprintf("%s_weekly_coordinator_action_list_%s.csv", org_name, current_date))
  excel_workbook_path <- file.path(output_dir, sprintf("%s_weekly_scored_appointments_%s.xlsx", org_name, current_date))
  
  # B. Style setup
  style_title     <- createStyle(fontName = "Calibri", fontSize = 16, fontColour = "#2F5496", textDecoration = "bold")
  style_subtitle  <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#595959", textDecoration = "italic")
  style_section   <- createStyle(fontName = "Calibri", fontSize = 12, fontColour = "#FFFFFF", fgFill = "#4472C4", textDecoration = "bold")
  style_header    <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#FFFFFF", fgFill = "#2F5496", textDecoration = "bold", halign = "center", valign = "center", border = "bottom", borderStyle = "medium")
  style_kpi_label <- createStyle(fontName = "Calibri", fontSize = 10, fontColour = "#2c3e50", textDecoration = "bold", halign = "right")
  style_kpi_value <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#2F5496", fgFill = "#F2F2F2", textDecoration = "bold", halign = "center", border = c("top", "bottom", "left", "right"), borderStyle = "thin")
  
  style_text_left <- createStyle(fontName = "Calibri", fontSize = 11, halign = "left")
  style_num_right <- createStyle(fontName = "Calibri", fontSize = 11, halign = "right")
  style_prob_pct  <- createStyle(fontName = "Calibri", fontSize = 11, halign = "right", numFmt = "0.0%")
  
  style_active_input <- createStyle(fontName = "Calibri", fontSize = 11, fgFill = "#DCE6F1", halign = "center", border = c("top", "bottom", "left", "right"), borderStyle = "thin")
  style_control_row  <- createStyle(fontName = "Calibri", fontSize = 11, fgFill = "#F2F2F2", fontColour = "#7F7F7F", halign = "center")
  
  style_active_datetime <- createStyle(
    fontName = "Calibri", fontSize = 11, fgFill = "#DCE6F1", halign = "center", 
    border = c("top", "bottom", "left", "right"), borderStyle = "thin", 
    numFmt = "yyyy-mm-dd hh:mm" 
  )
  
  # C. Export raw datasets
  write.csv(manifest, file = manifest_csv_path, row.names = FALSE)
  
  coordinator_list <- manifest %>%
    filter(trial_arm == "intervention") %>%
    select(appointment_id, patient_id, appointment_date, appointment_time, gp_practice, dna_risk) %>%
    arrange(desc(dna_risk))
  
  write.csv(coordinator_list, file = action_list_path, row.names = FALSE)
  
  # D. Generate styled Workbook
  message(sprintf("[%s] Generating clinical coordinator validation spreadsheet...", org_name))
  wb <- createWorkbook()
  
  # .........................................................................
  # SHEET 1: INSTRUCTIONS & MAPPING
  # .........................................................................
  sheet_name_inst <- "Instructions & mapping"
  addWorksheet(wb, sheet_name_inst)
  showGridLines(wb, sheet_name_inst, show = FALSE)
  
  writeData(wb, sheet_name_inst, "Weekly coordination export data specification", startCol = 1, startRow = 1)
  addStyle(wb, sheet_name_inst, style_title, rows = 1, cols = 1)
  writeData(wb, sheet_name_inst, "Operational guidelines and evaluation dictionary for clinical coordinators (v8)", startCol = 1, startRow = 2)
  addStyle(wb, sheet_name_inst, style_subtitle, rows = 2, cols = 1)
  
  writeData(wb, sheet_name_inst, "Model-generated variables vs manual logs", startCol = 1, startRow = 4)
  addStyle(wb, sheet_name_inst, style_section, rows = 4, cols = 1:3)
  setColWidths(wb, sheet_name_inst, cols = 1:3, widths = c(25, 50, 20))
  
  mapping_dict <- data.frame(
    Variable_Name = c(
      "appointment_id", "patient_id", "nhs_number", "appointment_date", "appointment_time", "appointment_day_of_week", "clinic_code", "tumoursite", "gp_practice",  
      "age", "sex", "ethnicity", "imd", "accessibility_flags", "dna_risk", "risk_label", "trial_arm",  
      "date_model_run", "model_version", "patient_landline_number", "patient_mobile_number",
      "tier1_text_timestamp", "tier2_text_timestamp", "t2_patient_intent", "tier3_phone_attempt", "tier3_call_timestamp",
      "number_of_call_attempts", "call_duration_minutes",
      "t3_patient_intent", "call_notes", "appointment_outcome"
    ),
    Description = c(
      "Appointment ID.", "Unique patient identifier.", "NHS number.", "Scheduled date of the appointment.", "Scheduled time of the appointment.", "Day of the week of the appointment.",
      "Clinic code.", "Primary tumour site/group.", "GP practice.", "Patient age.", "Patient sex.", "Ethnicity record.",
      "IMD decile.", "Active accessibility & vulnerability flags.", "Calibrated risk probability.",
      "Risk profile label.", "Randomised trial arm.", "Model execution date.", "Model version.",
      "Landline number.", "Mobile number.", "Tier 1 SMS timestamp.", "Tier 2 SMS timestamp.", "Patient intent derived from T2 interactive text.", "Phone outreach attempt.",
      "Phone call timestamp.", "Number of manual call attempts.", "Duration of the call (minutes).",
      "Patient intent derived from T3 coordinator call.", "Coordinator notes.", "Appointment outcome."
    ),
    Data_Source = rep("System/Log", 31)
  )
  
  writeData(wb, sheet_name_inst, t(c("Variable name", "Variable description", "Data source")), startCol = 1, startRow = 5, colNames = FALSE)
  addStyle(wb, sheet_name_inst, style_header, rows = 5, cols = 1:3)
  writeData(wb, sheet_name_inst, mapping_dict, startCol = 1, startRow = 6, colNames = FALSE)
  
  for (i in 6:(6 + nrow(mapping_dict) - 1)) {
    addStyle(wb, sheet_name_inst, style_text_left, rows = i, cols = 1:2)
    addStyle(wb, sheet_name_inst, style_text_left, rows = i, cols = 3)
  }
  
  # .........................................................................
  # SHEET 2: WEEKLY OUTREACH ROSTER (31 columns total)
  # .........................................................................
  sheet_name_rost <- "Weekly outreach roster"
  addWorksheet(wb, sheet_name_rost)
  showGridLines(wb, sheet_name_rost, show = FALSE)
  
  start_row <- 9
  end_row   <- start_row + nrow(manifest) - 1
  
  writeData(wb, sheet_name_rost, sprintf("Weekly coordination outreach roster: %s", org_name), startCol = 1, startRow = 1)
  addStyle(wb, sheet_name_rost, style_title, rows = 1, cols = 1)
  
  # KPI blocks referencing correct columns (risk_label is Column 16 / P, trial_arm is Column 17 / Q)
  writeData(wb, sheet_name_rost, "Total high-risk cohort:", startCol = 1, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 3, cols = 1)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(P%d:P%d, \"high-risk\")", start_row, end_row), startCol = 2, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 3, cols = 2)
  
  writeData(wb, sheet_name_rost, "Intervention arm (active calls):", startCol = 1, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 4, cols = 1)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(Q%d:Q%d, \"intervention\")", start_row, end_row), startCol = 2, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 4, cols = 2)
  
  writeData(wb, sheet_name_rost, "Control arm (passive standard care):", startCol = 4, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 3, cols = 4)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(Q%d:Q%d, \"control\")", start_row, end_row), startCol = 5, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 3, cols = 5)
  
  # Outreach completed is based on tier3_phone_attempt which is now column 25 / Y
  writeData(wb, sheet_name_rost, "Outreach completed:", startCol = 4, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 4, cols = 4)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(Y%d:Y%d, \"yes\")", start_row, end_row), startCol = 5, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 4, cols = 5)
  
  headers <- mapping_dict$Variable_Name
  
  writeData(wb, sheet_name_rost, t(headers), startCol = 1, startRow = 8, colNames = FALSE)
  addStyle(wb, sheet_name_rost, style_header, rows = 8, cols = 1:31)
  
  export_roster_data <- manifest %>%  
    select(
      appointment_id, patient_id, nhs_number, appointment_date, appointment_time, appointment_day_of_week, clinic_code, tumoursite, gp_practice,  
      age, sex, ethnicity, imd, accessibility_flags, dna_risk, risk_label, trial_arm,
      date_model_run, model_version, patient_landline_number, patient_mobile_number
    ) %>%
    mutate(
      tier1_text_timestamp = case_when(trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      tier2_text_timestamp = case_when(trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      t2_patient_intent = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      tier3_phone_attempt = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      tier3_call_timestamp = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      number_of_call_attempts = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      call_duration_minutes = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      t3_patient_intent = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      call_notes = case_when(trial_arm == "control" ~ "N/A - Control", trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ ""),
      appointment_outcome = case_when(trial_arm == "not in trial" ~ "N/A - Not in Trial", TRUE ~ "")
    )
  
  writeData(wb, sheet_name_rost, export_roster_data, startCol = 1, startRow = 9, colNames = FALSE)
  
  intervention_rows <- which(manifest$trial_arm == "intervention") + start_row - 1
  control_rows      <- which(manifest$trial_arm == "control") + start_row - 1
  not_in_trial_rows <- which(manifest$trial_arm == "not in trial") + start_row - 1
  
  # Base formatting style covering pre-populated metadata (Stopping at Column 21)
  addStyle(wb, sheet_name_rost, style_text_left, rows = start_row:end_row, cols = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 17, 18, 19, 20, 21), gridExpand = TRUE)
  addStyle(wb, sheet_name_rost, style_num_right, rows = start_row:end_row, cols = c(10, 13), gridExpand = TRUE) # age, imd
  addStyle(wb, sheet_name_rost, style_prob_pct, rows = start_row:end_row, cols = 15, gridExpand = TRUE) # dna_risk
  
  # Group-Based Vectorized Conditional Formatting
  if (length(intervention_rows) > 0) {
    # Blue background fields for manual logging columns
    addStyle(wb, sheet_name_rost, style_active_input, rows = intervention_rows, cols = c(24, 25, 27, 28, 29, 30, 31), gridExpand = TRUE)
    # Datetime styled columns for timestamps
    addStyle(wb, sheet_name_rost, style_active_datetime, rows = intervention_rows, cols = c(22, 23, 26), gridExpand = TRUE)
  }
  
  if (length(control_rows) > 0) {
    # Controls only receive passive texts and final outcome log
    addStyle(wb, sheet_name_rost, style_active_datetime, rows = control_rows, cols = c(22, 23), gridExpand = TRUE)
    addStyle(wb, sheet_name_rost, style_active_input, rows = control_rows, cols = 31, gridExpand = TRUE)
    # Locked/Greyed out phone tracking columns
    addStyle(wb, sheet_name_rost, style_control_row, rows = control_rows, cols = 24:30, gridExpand = TRUE)
  }
  
  if (length(not_in_trial_rows) > 0) {
    # Grey out manual logging fields completely
    addStyle(wb, sheet_name_rost, style_control_row, rows = not_in_trial_rows, cols = 22:31, gridExpand = TRUE)
  }
  
  # Pre-defined Dropdowns and Input Validations
  dataValidation(wb, sheet_name_rost, col = 24, rows = start_row:end_row, type = "list", value = '"confirm,cancel,reschedule,no_response"')
  dataValidation(wb, sheet_name_rost, col = 25, rows = start_row:end_row, type = "list", value = '"yes,no"')
  dataValidation(wb, sheet_name_rost, col = 27, rows = start_row:end_row, type = "whole", operator = "greaterThanOrEqual", value = "0")
  dataValidation(wb, sheet_name_rost, col = 28, rows = start_row:end_row, type = "whole", operator = "greaterThanOrEqual", value = "0")
  dataValidation(wb, sheet_name_rost, col = 29, rows = start_row:end_row, type = "list", value = '"confirm,cancel,reschedule,no_response"')
  dataValidation(wb, sheet_name_rost, col = 31, rows = start_row:end_row, type = "list", value = '"Attended,DNA,Cancelled,Rescheduled"')
  
  # Strict Decimal Datetime validation to block time-only "3pm" shortcuts
  # Removed custom validation (not supported in this openxlsx version)
  # Removed custom validation (not supported in this openxlsx version)
  # Removed custom validation (not supported in this openxlsx version)
  
  # Activate filters on all headers
  addFilter(wb, sheet_name_rost, row = 8, cols = 1:31)
  
  # Format Column Widths (Adjusted for Date and Time split columns)
  setColWidths(wb, sheet_name_rost, cols = 1, widths = 36)   # appointment_id
  setColWidths(wb, sheet_name_rost, cols = 2:3, widths = 14) # patient_id, nhs_number
  setColWidths(wb, sheet_name_rost, cols = 4, widths = 20)   # appointment_date (RENAMED)
  setColWidths(wb, sheet_name_rost, cols = 5, widths = 18)   # appointment_time (NEW)
  setColWidths(wb, sheet_name_rost, cols = 6, widths = 24)   # appointment_day_of_week
  setColWidths(wb, sheet_name_rost, cols = 7, widths = 14)   # clinic_code
  setColWidths(wb, sheet_name_rost, cols = 8, widths = 20)   # tumoursite
  setColWidths(wb, sheet_name_rost, cols = 9, widths = 26)   # gp_practice
  setColWidths(wb, sheet_name_rost, cols = 10, widths = 10)  # age
  setColWidths(wb, sheet_name_rost, cols = 11, widths = 10)  # sex
  setColWidths(wb, sheet_name_rost, cols = 12, widths = 16)  # ethnicity
  setColWidths(wb, sheet_name_rost, cols = 13, widths = 10)  # imd
  setColWidths(wb, sheet_name_rost, cols = 14, widths = 38)  # accessibility_flags
  setColWidths(wb, sheet_name_rost, cols = 15, widths = 14)  # dna_risk
  setColWidths(wb, sheet_name_rost, cols = 16, widths = 14)  # risk_label
  setColWidths(wb, sheet_name_rost, cols = 17, widths = 14)  # trial_arm
  setColWidths(wb, sheet_name_rost, cols = 18, widths = 22)  # date_model_run
  setColWidths(wb, sheet_name_rost, cols = 19, widths = 15)  # model_version
  setColWidths(wb, sheet_name_rost, cols = 20:21, widths = 20) # manual phone numbers (landline, mobile)
  setColWidths(wb, sheet_name_rost, cols = 22:23, widths = 26) # manual text timestamps (WIDENED)
  setColWidths(wb, sheet_name_rost, cols = 24, widths = 22)  # t2_patient_intent
  setColWidths(wb, sheet_name_rost, cols = 25:29, widths = 22) # manual call logging blanks (including call attempts, duration, t3 intent)
  setColWidths(wb, sheet_name_rost, cols = 30, widths = 30)  # call_notes
  setColWidths(wb, sheet_name_rost, cols = 31, widths = 22)  # appointment_outcome
  
  freezePane(wb, sheet_name_rost, firstActiveRow = 9, firstActiveCol = 2)
  dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
  saveWorkbook(wb, file = excel_workbook_path, overwrite = TRUE)
  
  # Console report
  n_total        <- nrow(manifest)
  n_high_risk    <- sum(manifest$risk_label == "high-risk")
  n_low_risk     <- sum(manifest$risk_label == "low-risk")
  n_control      <- sum(manifest$trial_arm == "control")
  n_intervention <- sum(manifest$trial_arm == "intervention")
  
  cat("\n====================================================\n")
  cat(sprintf("    Weekly prediction & allocation: %s        \n", org_name))
  cat("====================================================\n")
  cat(sprintf("Processed appointments:     %d patients\n", n_total))
  cat(sprintf("Operational threshold used: >= %.4f\n", op_threshold))
  cat("----------------------------------------------------\n")
  cat(sprintf("Low-risk (standard track):  %d patients (%.1f%%)\n", n_low_risk, (n_low_risk / n_total) * 100))
  cat(sprintf("High-risk (trial cohort):   %d patients (%.1f%%)\n", n_high_risk, (n_high_risk / n_total) * 100))
  cat("----------------------------------------------------\n")
  cat(sprintf("-> Control arm (passive):   %d patients\n", n_control))
  cat(sprintf("-> Intervention arm:        %d patients\n", n_intervention))
  cat("====================================================\n")
  cat(sprintf("Success. Scored manifest saved to %s\n", manifest_csv_path))
  cat(sprintf("Success. Pre-formatted coordinator workbook saved to %s\n", excel_workbook_path))
  cat("====================================================\n\n")
}


