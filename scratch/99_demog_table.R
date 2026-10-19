library(gtsummary)
library(cardx)
library(dplyr)
library(forcats)
library(flextable)
library(tidymodels)
library(yardstick)
library(probably)
library(broom)
library(config)

conf <- config::get()

# DATA SECTION
model_data_path <- "data/processed/train_engineered_500k_full.rds"
model_data <- readRDS(model_data_path)

# ==============================================================================
# 1. BACK-TRANSFORM & LUMP EXECUTIVE FEATURES
# ==============================================================================
table1_exec_data <- model_data %>%
  mutate(
    dna_outcome = fct_recode(dna_outcome, "Attended" = "attended"),
    fct_recode(dna_outcome, "attended" = "Attended"), 
    # A. Back-transform log lead time to real days (expm1 is inverse of log1p)
    lead_time_days = expm1(lead_time_days_log),
    
    # B. Lump granular age_group factor levels into 3 broad clinical stages
    age_broad = case_when(
      age_group %in% c("<18", "18-25", "18-24") ~ "< 25 years",
      age_group %in% c("25-34", "35-44", "45-54", "55-64") ~ "25–64 years",
      age_group %in% c("65-74", "75-84", "85-94", "95+") ~ "65+ years",
      TRUE ~ "Other / Unknown"
    ),
    age_broad = factor(age_broad, levels = c("< 25 years", "25–64 years", "65+ years", "Other / Unknown")),
    
    # C. Group 10 IMD Deciles into 5 Deprivation Quintiles
    imd_clean = as.character(imd),
    imd_quintile = case_when(
      imd_clean %in% c("1", "2")  ~ "Q1 (Most Deprived)",
      imd_clean %in% c("3", "4")  ~ "Q2",
      imd_clean %in% c("5", "6")  ~ "Q3",
      imd_clean %in% c("7", "8")  ~ "Q4",
      imd_clean %in% c("9", "10") ~ "Q5 (Least Deprived)",
      TRUE                        ~ "Unknown"
    ),
    imd_quintile = factor(imd_quintile, levels = c(
      "Q1 (Most Deprived)", "Q2", "Q3", "Q4", "Q5 (Least Deprived)", "Unknown"
    )),
    
    # D. Standardise Gender categories
    gender_clean = case_when(
      gender == "Female" ~ "Female",
      gender == "Male"   ~ "Male",
      TRUE               ~ "Other / Not Stated"
    ),
    
    # E. Simplify Appointment Type
    appt_type_clean = case_when(
      grepl("NEW", appointment_type, ignore.case = TRUE) ~ "New appointment",
      TRUE                                               ~ "Follow-up appointment"
    ),
    
    # F. Collapse all 20 'a_' vulnerability/accessibility flags into a single binary flag
    has_vulnerability = if_else(
      rowSums(across(starts_with("a_"), ~ .x == 1), na.rm = TRUE) > 0, 
      "Yes", "No"
    ),
    
    # G. Recode binary flags for clean table output
    is_morning_clean = if_else(is_morning == 1, "Yes", "No"),
    has_dna_history_clean = if_else(has_dna_history == 1, "Yes", "No")
  ) %>%
  # H. Lump high-cardinality predictors (GP practice, clinic code, specialties) to Top 5 + "Other"
  mutate(across(
    any_of(c("registered_gp_practice", "clinic_code", "site_code", "local_spec_code", "national_spec_code")),
    ~ fct_lump_n(factor(.), n = 3, other_level = "Other")
  ))

# ==============================================================================
# 2. DEFINE EXECUTIVE LABELS
# ==============================================================================
table1_labels <- list(
  age_broad              ~ "Age group",
  gender_clean           ~ "Sex / Gender",
  ethnicity_group        ~ "Ethnicity group",
  imd_quintile           ~ "Deprivation quintile (IMD)",
  distance_km            ~ "Travel distance (km)",
  lead_time_days         ~ "Appointment lead time (days)",
  is_morning_clean       ~ "Morning appointment slot",
  referral_urgency       ~ "Referral urgency",
  appt_type_clean        ~ "Appointment type",
  has_vulnerability      ~ "Has accessibility / vulnerability flag",
  registered_gp_practice ~ "Registered GP practice (top 3)",
  clinic_code            ~ "Clinic code (top 3)",
  national_spec_code     ~ "Speciality (top 3)",
  has_dna_history_clean  ~ "Prior DNA history (past 12 months)"
)

# ==============================================================================
# 3. BUILD CONCISE 1-PAGE SUMMARY TABLE
# ==============================================================================
tbl_executive <- table1_exec_data %>%
  select(
    dna_outcome, age_broad, gender_clean, ethnicity_group, 
    imd_quintile, distance_km, lead_time_days, is_morning_clean, 
    referral_urgency, appt_type_clean, has_vulnerability, 
    registered_gp_practice, clinic_code, national_spec_code, has_dna_history_clean
  ) %>%
  tbl_summary(
    by = dna_outcome,
    label = table1_labels,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      lead_time_days   ~ 0 # Round lead time days to whole numbers
    ),
    missing_text = "Unknown"
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(test = list(
    all_categorical() ~ "chisq.test",
    all_continuous()  ~ "wilcox.test"
  )) %>%
  bold_labels()



# ==============================================================================
# 4. EXPORT DIRECTLY TO WORD
# ==============================================================================
tbl_executive %>%
  as_flex_table() %>%
  autofit() %>%
  save_as_docx(path = "outputs/table1_executive_summary.docx")




#MODEL SECTION

# ==============================================================================
# MODEL METRICS RE-EVALUATION SUITE
# ==============================================================================
# Evaluates discrimination, calibration, statistical rank separation, 
# and operational lift directly from `model_bundle` and dataset.



model_bundle_path <- "data/processed/models/rf_final_model_500k_full.rds"
model_bundle  <- readRDS(model_bundle_path)


evaluate_model_bundle <- function(
    model_bundle_path = "data/processed/models/rf_final_model_500k_full.rds",
    op_threshold      = 0.10
) {
  
  message("--- Loading Model Bundle ---")
  if (!file.exists(model_bundle_path)) stop("Model bundle file not found: ", model_bundle_path)
  
  model_bundle <- readRDS(model_bundle_path)
  
  # 0. PR/ROC plots
  
  # -------------------------------------------------------------------------
  # C. Generate ROC & PR Plots (Using the safe best_preds) [cite: 6]
  # -------------------------------------------------------------------------
  roc_plot <- model_bundle$predictions |> (\(preds) {
    auc_val <- preds |> roc_auc(truth = !!sym(conf$target_col), .pred_DNA) |> pull(.estimate)
    preds |>
      roc_curve(truth = !!sym(conf$target_col), .pred_DNA) |>
      ggplot(aes(x = 1 - specificity, y = sensitivity)) +
      geom_abline(slope = 1, linetype = 2, colour = "red") +
      geom_path(linewidth = 1, color = "midnightblue") +
      coord_equal() +
      theme_bw() +
      labs(title = "ROC curve", subtitle = paste0("ROC AUC: ", round(auc_val, 3)))
  })()
  
  pr_plot <- model_bundle$predictions |> (\(preds) {
    pr_val <- preds |> pr_auc(truth = !!sym(conf$target_col), .pred_DNA) |> pull(.estimate)
    preds |>
      pr_curve(truth = !!sym(conf$target_col), .pred_DNA) %>%
      ggplot(aes(x = recall, y = precision)) +
      geom_path(linewidth = 1, color = "midnightblue") +
      geom_hline(yintercept = 0.05, lty = 2, color = "red") +
      coord_equal() +
      theme_bw() +
      labs(title = "PR Curve", subtitle = paste0("PR AUC: ", round(pr_val, 3)))
  })()
  
  pr_roc_ptc <- patchwork::wrap_plots(list(pr_plot, roc_plot));pr_roc_ptc
  
  # 1. Harvest raw out-of-fold CV predictions
  cv_preds   <- model_bundle$predictions
  calibrator <- model_bundle$calibrator
  
  # 2. Preserve raw probabilities BEFORE applying calibrator
  cv_preds <- cv_preds %>%
    mutate(.pred_DNA_raw = .pred_DNA) %>%
    cal_apply(calibrator) # Overwrites .pred_DNA with calibrated values
  
  # 3. Calculate Cohort Totals & Baseline Prevalence
  n_total <- nrow(cv_preds)
  n_dna   <- sum(cv_preds$dna_outcome == "DNA", na.rm = TRUE)
  p_base  <- n_dna / n_total
  
  # 4. Discrimination & Calibration Metrics
  roc_auc_val <- cv_preds %>% 
    roc_auc(truth = dna_outcome, .pred_DNA, event_level = "first") %>% 
    pull(.estimate)
  
  pr_auc_val <- cv_preds %>% 
    pr_auc(truth = dna_outcome, .pred_DNA, event_level = "first") %>% 
    pull(.estimate)
  
  brier_score <- mean((if_else(cv_preds$dna_outcome == "DNA", 1, 0) - cv_preds$.pred_DNA)^2, na.rm = TRUE)
  
  # 5. Mann-Whitney U Test & Score Summaries
  wilcox_res <- wilcox.test(.pred_DNA ~ dna_outcome, data = cv_preds, exact = FALSE)
  
  score_summary <- cv_preds %>%
    group_by(dna_outcome) %>%
    summarise(
      median_p = median(.pred_DNA, na.rm = TRUE),
      p25      = quantile(.pred_DNA, 0.25, na.rm = TRUE),
      p75      = quantile(.pred_DNA, 0.75, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # 6. Top 10% Decile Lift
  decile_cutoff      <- quantile(cv_preds$.pred_DNA, 0.90, na.rm = TRUE)
  top_decile_preds   <- cv_preds %>% filter(.pred_DNA >= decile_cutoff)
  top_decile_tp      <- sum(top_decile_preds$dna_outcome == "DNA", na.rm = TRUE)
  top_decile_prec    <- top_decile_tp / nrow(top_decile_preds)
  top_decile_lift    <- top_decile_prec / p_base
  top_decile_capture <- top_decile_tp / n_dna
  
  # 7. Operational Threshold Lift (At op_threshold)
  flagged_preds <- cv_preds %>% filter(.pred_DNA >= op_threshold)
  n_flagged     <- nrow(flagged_preds)
  op_tp         <- sum(flagged_preds$dna_outcome == "DNA", na.rm = TRUE)
  op_prec       <- if_else(n_flagged > 0, op_tp / n_flagged, 0)
  op_lift       <- op_prec / p_base
  op_capture    <- op_tp / n_dna
  
  # 8. Print Executive Console Summary
  cat("\n==================================================\n")
  cat("               EVALUATION METRICS SUMMARY REPORT                   \n")
  cat("==================================================\n")
  cat(sprintf("Cohort Total Volume (N):       %s appointments\n", format(n_total, big.mark = ",")))
  cat(sprintf("Observed Total DNAs:           %s (%0.2f%% baseline rate)\n", format(n_dna, big.mark = ","), p_base * 100))
  cat("--------------------------------------------------\n")
  cat("1. MODEL DISCRIMINATION & CALIBRATION:\n")
  cat(sprintf("   - Area Under PR (PR-AUC):              %0.4f (Random Baseline: %0.4f)\n", pr_auc_val, p_base))
  cat(sprintf("   - Area Under ROC (ROC-AUC):             %0.4f\n", roc_auc_val))
  cat(sprintf("   - Brier Score (Post-Platt):             %0.4f\n", brier_score))
  cat("--------------------------------------------------\n")
  cat("2. RANK SEPARATION (MANN-WHITNEY U TEST):\n")
  cat(sprintf("   - Wilcoxon p-value:                     %s\n", format.pval(wilcox_res$p.value, digits = 3)))
  cat(sprintf("   - Common Language Effect Size (CLES):   %0.2f%%\n", roc_auc_val * 100))
  for (i in 1:nrow(score_summary)) {
    cat(sprintf("   - %-10s Risk Median (IQR):      %0.3f (%0.3f - %0.3f)\n", 
                score_summary$dna_outcome[i], score_summary$median_p[i], score_summary$p25[i], score_summary$p75[i]))
  }
  cat("--------------------------------------------------\n")
  cat(sprintf("3. OPERATIONAL THRESHOLD (>= %0.2f):\n", op_threshold))
  cat(sprintf("   - Population Flagged:                   %s (%0.2f%% of clinic)\n", format(n_flagged, big.mark = ","), (n_flagged / n_total) * 100))
  cat(sprintf("   - Precision (Hit Rate):                 %0.2f%%\n", op_prec * 100))
  cat(sprintf("   - Operational Lift:                     %0.2fx\n", op_lift))
  cat(sprintf("   - Total DNA Capture Rate (Recall):      %0.2f%%\n", op_capture * 100))
  cat("--------------------------------------------------\n")
  cat("4. TOP 10% DECILE PERFORMANCE:\n")
  cat(sprintf("   - Top Decile Score Cutoff:              >= %0.4f\n", decile_cutoff))
  cat(sprintf("   - Top Decile Precision:                 %0.2f%%\n", top_decile_prec * 100))
  cat(sprintf("   - Top Decile Operational Lift:          %0.2fx\n", top_decile_lift))
  cat(sprintf("   - Top Decile DNA Capture Rate:          %0.2f%%\n", top_decile_capture * 100))
  cat("==================================================\n\n")
  
  return(invisible(list(
    n_total = n_total, n_dna = n_dna, p_base = p_base,
    pr_auc = pr_auc_val, roc_auc = roc_auc_val, brier_score = brier_score,
    top_decile_prec = top_decile_prec, top_decile_lift = top_decile_lift,
    cv_preds = cv_preds, pr_roc_ptc = pr_roc_ptc
  )))
}



# Execute evaluation
metrics_summary <- evaluate_model_bundle(
  model_bundle_path = "data/processed/models/rf_final_model_500k_full.rds",
  op_threshold      = 0.10
)




ggsave(filename = "materials/pr_roc.png",
       metrics_summary$pr_roc_ptc,
       width = 10,
       height = 6,
       scale = 0.6,
       bg = "white")
