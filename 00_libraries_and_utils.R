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

apply_custom_feature_engineering <- function(data) {
  # Handle IMD coalescing from the raw database column if present
  if ("index_multiple_deprivation_decile" %in% names(data)) {
    data <- data %>%
      mutate(
        imd = coalesce(
          as.character(index_multiple_deprivation_decile),
          "unknown"
        )
      )
  }
  
  data %>%
    mutate(
      ethnicity_clean = collect_ethnicity(ethnicity),
      ethnicity_group = aggregate_ethnicity_high_level(ethnicity_clean),
      # Set the baseline reference level to the most frequent category
      ethnicity_group = factor(ethnicity_group) %>% relevel(ref = "White"),
      appt_date = as.Date(substring(appt_month, 1, 10), format = "%d/%m/%Y"),
      appt_dow = factor(weekdays(appt_date)),
      appt_month_num = as.factor(format(appt_date, "%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      lead_time_days_log = log1p(pmax(0, lead_time_days)),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      appt_hour_sin = sin(2 * pi * appt_hour / 24),
      appt_hour_cos = cos(2 * pi * appt_hour / 24),
      has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
    ) %>%
    # Drop raw redundant columns safely
    select(
      -any_of(c(
        "appt_hour", "lead_time_days", "appt_date", "appt_month", 
        "prev_dna_ly", "ethnicity", "ethnicity_clean", "age_group",
        "index_multiple_deprivation_decile"
      ))
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


# =========================================================================
# OPERATIONAL EXPORT: EXCEL WRITER WITH VALIDATION
# =========================================================================
# This function takes the scored weekly manifest and writes it directly
# into a beautifully styled Excel workbook with dynamic data validations.

write_weekly_coordination_sheet <- function(final_manifest, output_path = "outputs/weekly_scored_appointments.xlsx") {
  
  # Ensure the openxlsx package is available
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    install.packages("openxlsx")
  }
  library(openxlsx)
  library(dplyr)
  
  # Create a fresh workbook
  wb <- createWorkbook()
  
  # Define professional style sheets (navy, light blue, and light gray)
  style_title <- createStyle(fontName = "Calibri", fontSize = 16, fontColour = "#2F5496", textDecoration = "bold")
  style_subtitle <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#595959", italic = TRUE)
  style_section <- createStyle(fontName = "Calibri", fontSize = 12, fontColour = "#FFFFFF", fgFill = "#4472C4", textDecoration = "bold")
  style_header <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#FFFFFF", fgFill = "#2F5496", textDecoration = "bold", halign = "center", valign = "center", border = "bottom", borderStyle = "medium")
  
  style_kpi_label <- createStyle(fontName = "Calibri", fontSize = 10, fontColour = "#2c3e50", textDecoration = "bold", halign = "right")
  style_kpi_value <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#2F5496", fgFill = "#F2F2F2", textDecoration = "bold", halign = "center", border = "surrounding", borderStyle = "thin")
  
  style_text_left <- createStyle(fontName = "Calibri", fontSize = 11, halign = "left")
  style_num_right <- createStyle(fontName = "Calibri", fontSize = 11, halign = "right")
  style_prob_pct  <- createStyle(fontName = "Calibri", fontSize = 11, halign = "right", numFmt = "0.0%")
  
  style_active_input <- createStyle(fontName = "Calibri", fontSize = 11, fgFill = "#DCE6F1", halign = "center", border = "surrounding", borderStyle = "thin")
  style_control_row  <- createStyle(fontName = "Calibri", fontSize = 11, fgFill = "#F2F2F2", fontColour = "#7F7F7F", halign = "center")
  
  # -------------------------------------------------------------------------
  # SHEET 1: INSTRUCTIONS & MAPPING
  # -------------------------------------------------------------------------
  sheet_name_inst <- "Instructions & mapping"
  addWorksheet(wb, sheet_name_inst)
  showGridLines(wb, sheet_name_inst, show = FALSE)
  
  # Add titles
  writeData(wb, sheet_name_inst, "Weekly coordination export data specification", startCol = 1, startRow = 1)
  addStyle(wb, sheet_name_inst, style_title, rows = 1, cols = 1)
  writeData(wb, sheet_name_inst, "Operational guidelines and evaluation dictionary for clinical coordinators", startCol = 1, startRow = 2)
  addStyle(wb, sheet_name_inst, style_subtitle, rows = 2, cols = 1)
  
  # Section: Variables mapping
  writeData(wb, sheet_name_inst, "Model-generated variables vs manual logs", startCol = 1, startRow = 4)
  addStyle(wb, sheet_name_inst, style_section, rows = 4, cols = 1:3)
  setColWidths(wb, sheet_name_inst, cols = 1:3, widths = c(25, 45, 20))
  
  # Build mapping dictionary table
  mapping_dict <- data.frame(
    Variable_Name = c(
      "patient_id", "age_at_appointment", "registered_gp_practice", "imd", 
      "appt_dow", "dna_probability", "trial_arm", "phone_contact_occurred", 
      "call_timestamp", "patient_intent", "cancel_reschedule_wish", "cancellation_reason"
    ),
    Description = c(
      "Unique patient identifier generated by system.",
      "Patient age on the scheduled appointment date.",
      "GP practice where the patient is currently registered.",
      "Index of Multiple Deprivation decile (1-10 or unknown).",
      "Calculated day of the week for the appointment.",
      "Calibrated risk probability score of a no-show.",
      "Randomized trial allocation arm.",
      "Manual log: did a successful telephone contact occur?",
      "Manual log: precise date and time of the outreach attempt.",
      "Manual log: patient's stated intent (confirm, cancel, reschedule).",
      "Manual log: did the patient express a wish to cancel or reschedule?",
      "Manual log: primary reason for cancellation if applicable."
    ),
    Data_Source = c(
      "Model generated", "EHR record", "EHR record", "Census lookup",
      "System computed", "Model generated", "Model generated", "Coordinator log",
      "Coordinator log", "Coordinator log", "Coordinator log", "Coordinator log"
    )
  )
  
  # Write table headers and content
  writeData(wb, sheet_name_inst, t(c("Variable name", "Variable description", "Data source")), startCol = 1, startRow = 5, colNames = FALSE)
  addStyle(wb, sheet_name_inst, style_header, rows = 5, cols = 1:3)
  writeData(wb, sheet_name_inst, mapping_dict, startCol = 1, startRow = 6, colNames = FALSE)
  for (i in 6:(6 + nrow(mapping_dict) - 1)) {
    addStyle(wb, sheet_name_inst, style_text_left, rows = i, cols = 1:2)
    addStyle(wb, sheet_name_inst, style_text_left, rows = i, cols = 3)
  }
  
  # Section: Downstream evaluation map
  row_eval <- 20
  writeData(wb, sheet_name_inst, "Evaluation model mapping schema", startCol = 1, startRow = row_eval)
  addStyle(wb, sheet_name_inst, style_section, rows = row_eval, cols = 1:3)
  
  eval_notes <- c(
    "Model 1: Intent-to-Treat (ITT) Advance Cancellations",
    "• Evaluated on all randomized patients (including DNAs and attendances).",
    "• Evaluates the global impact of the text intervention program on freeing up clinic slots.",
    "• Core metrics: 'trial_arm' predicts 'cancellation_outcome'.",
    "",
    "Model 2: Per-Protocol (PP) Wasted Capacity / DNA",
    "• Evaluated strictly on at-risk patients (excludes successful advance cancellations).",
    "• Evaluates the specific behavioral impact of automated texts vs. human phone outreach.",
    "• Core metrics: 'phone_contact_occurred' segments exposure tiers; 'dna_probability' is used as a log-offset."
  )
  for (i in seq_along(eval_notes)) {
    writeData(wb, sheet_name_inst, eval_notes[i], startCol = 1, startRow = row_eval + i)
    addStyle(wb, sheet_name_inst, style_text_left, rows = row_eval + i, cols = 1)
  }
  
  # Section: Data integrity guidelines
  row_integrity <- 31
  writeData(wb, sheet_name_inst, "Clinical coordinator data integrity guidelines", startCol = 1, startRow = row_integrity)
  addStyle(wb, sheet_name_inst, style_section, rows = row_integrity, cols = 1:3)
  
  integrity_notes <- c(
    "1. Telephony logging synchronization:",
    "• Always verify the patient ID in the VoIP system matches the spreadsheet before making calls.",
    "2. Avoid cancellation misclassifications:",
    "• Do not delete appointments or mark them as DNA if the patient verbally cancels during outreach.",
    "3. Time sequence accuracy:",
    "• Log the call timestamp immediately. Back-dating call records corrupts the latency evaluation."
  )
  for (i in seq_along(integrity_notes)) {
    writeData(wb, sheet_name_inst, integrity_notes[i], startCol = 1, startRow = row_integrity + i)
    addStyle(wb, sheet_name_inst, style_text_left, rows = row_integrity + i, cols = 1)
  }
  
  # -------------------------------------------------------------------------
  # SHEET 2: WEEKLY OUTREACH ROSTER
  # -------------------------------------------------------------------------
  sheet_name_rost <- "Weekly outreach roster"
  addWorksheet(wb, sheet_name_rost)
  showGridLines(wb, sheet_name_rost, show = FALSE)
  
  # Dynamic indexing of rows
  start_row <- 9
  end_row   <- start_row + nrow(final_manifest) - 1
  
  # Add title block
  writeData(wb, sheet_name_rost, "Weekly coordination outreach roster", startCol = 1, startRow = 1)
  addStyle(wb, sheet_name_rost, style_title, rows = 1, cols = 1)
  
  # Add summary cards (KPI Block)
  writeData(wb, sheet_name_rost, "Total high-risk cohort:", startCol = 1, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 3, cols = 1)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTA(A%d:A%d)", start_row, end_row), startCol = 2, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 3, cols = 2)
  
  writeData(wb, sheet_name_rost, "Intervention arm (active calls):", startCol = 1, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 4, cols = 1)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(G%d:G%d, \"intervention\")", start_row, end_row), startCol = 2, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 4, cols = 2)
  
  writeData(wb, sheet_name_rost, "Control arm (passive standard care):", startCol = 4, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 3, cols = 4)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(G%d:G%d, \"control\")", start_row, end_row), startCol = 5, startRow = 3)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 3, cols = 5)
  
  writeData(wb, sheet_name_rost, "Outreach completed:", startCol = 4, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_label, rows = 4, cols = 4)
  writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(H%d:H%d, \"yes\")", start_row, end_row), startCol = 5, startRow = 4)
  addStyle(wb, sheet_name_rost, style_kpi_value, rows = 4, cols = 5)
  
  # Column headers
  headers <- c(
    "patient_id", "age_at_appointment", "registered_gp_practice", "imd", 
    "appt_dow", "dna_probability", "trial_arm", "phone_contact_occurred", 
    "call_timestamp", "patient_intent", "cancel_reschedule_wish", "cancellation_reason"
  )
  writeData(wb, sheet_name_rost, t(headers), startCol = 1, startRow = 8, colNames = FALSE)
  addStyle(wb, sheet_name_rost, style_header, rows = 8, cols = 1:12)
  
  # Populating the scored model output variables
  writeData(wb, sheet_name_rost, final_manifest %>% select(patient_id, age_at_appointment, registered_gp_practice, imd, appt_dow, dna_probability, trial_arm), startCol = 1, startRow = 9, colNames = FALSE)
  
  # Apply formatting to data rows and populate coordinator fields dynamically
  for (r in start_row:end_row) {
    # Extract the trial arm for the current row
    arm_value <- final_manifest$trial_arm[r - start_row + 1]
    
    # 1. Base model predictions formatting (Columns A-G)
    addStyle(wb, sheet_name_rost, style_text_left, rows = r, cols = c(1, 3, 7))
    addStyle(wb, sheet_name_rost, style_num_right, rows = r, cols = c(2, 4))
    addStyle(wb, sheet_name_rost, style_text_left, rows = r, cols = 5)
    addStyle(wb, sheet_name_rost, style_prob_pct, rows = r, cols = 6)
    
    # 2. Coordinator interactive logging formatting (Columns H-L)
    if (arm_value == "intervention") {
      # Active Intervention patients: Highlight fields in light blue ("edit me" signal)
      addStyle(wb, sheet_name_rost, style_active_input, rows = r, cols = 8:12)
    } else {
      # Control group: pre-populate with "N/A - Control" and gray out to prevent edits
      writeData(wb, sheet_name_rost, "N/A - Control", startCol = 8, startRow = r)
      writeData(wb, sheet_name_rost, "N/A - Control", startCol = 9, startRow = r)
      writeData(wb, sheet_name_rost, "N/A - Control", startCol = 10, startRow = r)
      writeData(wb, sheet_name_rost, "N/A - Control", startCol = 11, startRow = r)
      writeData(wb, sheet_name_rost, "N/A - Control", startCol = 12, startRow = r)
      addStyle(wb, sheet_name_rost, style_control_row, rows = r, cols = 8:12)
    }
  }
  
  # Apply dynamic data validation dropdown lists for active intervention rows
  dataValidation(wb, sheet_name_rost, col = 8, rows = start_row:end_row, type = "list", value = '\"yes,no\"')
  dataValidation(wb, sheet_name_rost, col = 10, rows = start_row:end_row, type = "list", value = '\"confirm,cancel,reschedule,no_response\"')
  dataValidation(wb, sheet_name_rost, col = 11, rows = start_row:end_row, type = "list", value = '\"yes,no\"')
  dataValidation(wb, sheet_name_rost, col = 12, rows = start_row:end_row, type = "list", value = '\"Scheduling conflict,Forgot or misplaced details,Transport or mobility barrier,Unwell or acute clinical issue,Appointment no longer required,Declined to state / other\"')
  
  # Column width padding
  setColWidths(wb, sheet_name_rost, cols = 1, widths = 14)   # patient_id
  setColWidths(wb, sheet_name_rost, cols = 2, widths = 20)   # age_at_appointment
  setColWidths(wb, sheet_name_rost, cols = 3, widths = 26)   # gp_practice
  setColWidths(wb, sheet_name_rost, cols = 4, widths = 10)   # imd
  setColWidths(wb, sheet_name_rost, cols = 5, widths = 14)   # appt_dow
  setColWidths(wb, sheet_name_rost, cols = 6, widths = 16)   # dna_probability
  setColWidths(wb, sheet_name_rost, cols = 7, widths = 14)   # trial_arm
  setColWidths(wb, sheet_name_rost, cols = 8:12, widths = 22) # manual blanks
  
  # Lock split panes below header row so scrolling is coordinator-friendly
  freezePanes(wb, sheet_name_rost, firstActiveRow = 9, firstActiveCol = 1)
  
  # Save workbook
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  saveWorkbook(wb, file = output_path, overwrite = TRUE)
}
