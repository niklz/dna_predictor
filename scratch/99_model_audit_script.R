# =========================================================================
# SCRIPT: MODEL PREDICTION AUDIT & DEBIASED VALIDATION
# =========================================================================
# This script performs a rigorous two-part audit to verify model stability:
# 
# 1. Model Drift Audit:
#    Compares prediction distribution, individual risk scores, and clinical
#    high-risk classifications between the legacy model and the new model.
#
# 2. Feature Ablation Audit (Month Number Impact):
#    Directly measures how much the date-parsing fix (which restored month_num)
#    altered predictions compared to leaving the column empty (as NA/Unknown).
# -------------------------------------------------------------------------

source("00_libraries_and_utils.R")
conf <- config::get()

# Set up output directory
dir.create("outputs/audit", showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Load Evaluation Data & Models
# -------------------------------------------------------------------------
message("Loading audit datasets...")
# We use the test/evaluation dataset to guarantee unbiased comparisons
test_data <- tryCatch({
  readRDS("data/processed/test_raw.rds")
}, error = function(e) {
  # Fallback to train engineered if test raw is not saved separately
  message("test_raw.rds not found. Using train_engineered.rds for audit...")
  readRDS("data/processed/train_engineered.rds")
})

# Load the new finalized model bundle
new_model_path <- "data/processed/rf_final_model.rds"
if (!file.exists(new_model_path)) {
  stop("New model bundle not found. Please run your training pipeline first.")
}
new_bundle <- readRDS(new_model_path)
new_model  <- new_bundle$model
new_calibrator <- new_bundle$calibrator

# Load the legacy model bundle
# NOTE: Ensure you copy your old model to this path to run the comparative audit!
legacy_model_path <- "data/processed/rf_final_model_legacy.rds"
has_legacy_model  <- file.exists(legacy_model_path)

if (has_legacy_model) {
  legacy_bundle <- readRDS(legacy_model_path)
  legacy_model  <- legacy_bundle$model
  legacy_calibrator <- legacy_bundle$calibrator
} else {
  message("ℹ Legacy model not found at data/processed/rf_final_model_legacy.rds")
  message("ℹ Skipping direct Model-vs-Model comparison. Running Month Ablation Audit only.")
}

# -------------------------------------------------------------------------
# 2. Part 1: Model-vs-Model Drift Audit (GLM vs Mixed, etc.)
# -------------------------------------------------------------------------
if (has_legacy_model) {
  message("\n--- Running Model-vs-Model Drift Audit ---")
  
  # A. Feature Engineer raw test data for both models (just in case recipe changed)
  test_prepped_new <- apply_custom_feature_engineering(test_data)
  test_prepped_legacy <- apply_custom_feature_engineering(test_data)
  
  # B. Generate Calibrated predictions from both models
  pred_new <- predict(new_model, new_data = test_prepped_new, type = "prob") %>%
    cal_apply(new_calibrator) %>%
    pull(.pred_DNA)
    
  pred_legacy <- predict(legacy_model, new_data = test_prepped_legacy, type = "prob") %>%
    cal_apply(legacy_calibrator) %>%
    pull(.pred_DNA)
    
  # C. Calculate Distribution Metrics
  mae_drift  <- mean(abs(pred_new - pred_legacy))
  rmse_drift <- sqrt(mean((pred_new - pred_legacy)^2))
  max_drift  <- max(abs(pred_new - pred_legacy))
  corr_pears <- cor(pred_new, pred_legacy, method = "pearson")
  corr_spear <- cor(pred_new, pred_legacy, method = "spearman")
  
  # D. Load operational threshold to check classification agreement
  threshold_data <- readRDS("data/processed/risk_threshold.RDS")
  op_threshold   <- threshold_data$production_threshold
  
  class_new    <- ifelse(pred_new >= op_threshold, "High-Risk", "Low-Risk")
  class_legacy <- ifelse(pred_legacy >= op_threshold, "High-Risk", "Low-Risk")
  
  classification_agreement <- mean(class_new == class_legacy)
  
  # Calculate high-risk cohort intersection (Operational Overlap)
  high_risk_new_idx <- which(class_new == "High-Risk")
  high_risk_leg_idx <- which(class_legacy == "High-Risk")
  
  intersection_count <- length(intersect(high_risk_new_idx, high_risk_leg_idx))
  union_count        <- length(union(high_risk_new_idx, high_risk_leg_idx))
  jaccard_overlap    <- intersection_count / union_count
  
  # Report Core Metrics to Console
  cat("\n====================================================\n")
  cat("          MODEL DRIFT AUDIT REPORT                  \n")
  cat("====================================================\n")
  cat(sprintf("Core Drift Statistics:\n"))
  cat(sprintf("- Mean Absolute Drift:       %.4f (Avg score change)\n", mae_drift))
  cat(sprintf("- Root Mean Squared Drift:   %.4f\n", rmse_drift))
  cat(sprintf("- Max Individual Drift:      %.4f\n", max_drift))
  cat(sprintf("- Pearson Correlation:       %.4f (Linear rank agreement)\n", corr_pears))
  cat(sprintf("- Spearman Rank Correlation: %.4f (Ordinal rank agreement)\n", corr_spear))
  cat("----------------------------------------------------\n")
  cat(sprintf("Operational Impact (Threshold >= %.4f):\n", op_threshold))
  cat(sprintf("- Global Class Agreement:    %.2f%%\n", classification_agreement * 100))
  cat(sprintf("- High-Risk Jaccard Overlap: %.2f%% (Shared coordinator pool)\n", jaccard_overlap * 100))
  cat(sprintf("  * Legacy model flagged:     %d patients\n", length(high_risk_leg_idx)))
  cat(sprintf("  * New model flagged:        %d patients\n", length(high_risk_new_idx)))
  cat(sprintf("  * Shared patients:          %d patients\n", intersection_count))
  cat("====================================================\n\n")
  
  # E. Generate Visual Diagnostics
  # 1. Density overlap plot
  compare_df <- tibble(
    score = c(pred_legacy, pred_new),
    model = rep(c("Legacy Model (Mixed Encoding)", "New Model (GLM Encoding)"), each = length(pred_new))
  )
  
  p_density <- ggplot(compare_df, aes(x = score, fill = model)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = c("#7f8c8d", "#2980b9")) +
    labs(
      title = "Risk Score Density Distribution Shift",
      subtitle = "Compares overall probability distributions between model versions",
      x = "Calibrated DNA Risk Probability",
      y = "Density",
      fill = "Model Version"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave("outputs/audit/model_density_shift.png", plot = p_density, width = 8, height = 5)
  
  # 2. Scatter Correlation plot
  p_scatter <- tibble(Legacy = pred_legacy, New = pred_new) %>%
    ggplot(aes(x = Legacy, y = New)) +
    geom_point(alpha = 0.1, color = "midnightblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = "Individual Risk Score Correlation",
      subtitle = sprintf("Pearson Correlation: %.4f | Red line represents perfect agreement", corr_pears),
      x = "Legacy Calibrated Score",
      y = "New Calibrated Score"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("outputs/audit/model_score_correlation.png", plot = p_scatter, width = 6, height = 6)
  
  message("Model comparative audit completed. Diagnostic charts saved to 'outputs/audit/'")
}

# -------------------------------------------------------------------------
# 3. Part 2: Feature Ablation Audit (Month Number Impact)
# -------------------------------------------------------------------------
message("\n--- Running Month Number Feature Ablation Audit ---")

# Step A: Create month-corrupted data mimicking old date-parsing bug
# (Forces appt_month_num and appt_dow to be NA as it was last week)
test_engineered_clean <- apply_custom_feature_engineering(test_data)

test_corrupted <- test_data %>%
  mutate(appt_month = NA_character_) # Wipes appt_month prior to feature engineering

test_engineered_corrupted <- apply_custom_feature_engineering(test_corrupted)

# Step B: Score both datasets using the *NEW* model
pred_clean <- predict(new_model, new_data = test_engineered_clean, type = "prob") %>%
  cal_apply(new_calibrator) %>%
  pull(.pred_DNA)

pred_corrupted <- predict(new_model, new_data = test_engineered_corrupted, type = "prob") %>%
  cal_apply(new_calibrator) %>%
  pull(.pred_DNA)

# Step C: Calculate the Delta introduced purely by Month/Date Features
ablation_mae  <- mean(abs(pred_clean - pred_corrupted))
ablation_max  <- max(abs(pred_clean - pred_corrupted))
ablation_corr <- cor(pred_clean, pred_corrupted, method = "spearman")

# Step D: Operational Threshold Shift Analysis
threshold_data <- readRDS("data/processed/risk_threshold.RDS")
op_threshold   <- threshold_data$production_threshold

class_clean     <- ifelse(pred_clean >= op_threshold, "High-Risk", "Low-Risk")
class_corrupted <- ifelse(pred_corrupted >= op_threshold, "High-Risk", "Low-Risk")

overlap_rate <- mean(class_clean == class_corrupted)

high_risk_clean_idx     <- which(class_clean == "High-Risk")
high_risk_corrupted_idx <- which(class_corrupted == "High-Risk")

intersection_ablation <- length(intersect(high_risk_clean_idx, high_risk_corrupted_idx))
union_ablation        <- length(union(high_risk_clean_idx, high_risk_corrupted_idx))
jaccard_ablation      <- intersection_ablation / union_ablation

# Report Month-Num Impact to Console
cat("====================================================\n")
cat("      MONTH FEATURE ABLATION AUDIT REPORT           \n")
cat("====================================================\n")
cat(sprintf("Predictive Impact of Date Parsing Bug:\n"))
cat(sprintf("- Mean Score Deviation (MAE): %.4f (Average absolute risk shift)\n", ablation_mae))
cat(sprintf("- Maximum Score Deviation:    %.4f (Largest shift for one patient)\n", ablation_max))
cat(sprintf("- Spearman Rank Correlation:  %.4f (Ordinal rank agreement)\n", ablation_corr))
cat("----------------------------------------------------\n")
cat(sprintf("Operational Impact on Trial Allocation (Threshold >= %.4f):\n", op_threshold))
cat(sprintf("- Roster Overlap Rate:        %.2f%%\n", overlap_rate * 100))
cat(sprintf("- High-Risk Pool Jaccard:     %.2f%% (Shared cohort agreement)\n", jaccard_ablation * 100))
cat(sprintf("  * Scored with month data:    %d patients\n", length(high_risk_clean_idx)))
cat(sprintf("  * Scored without month data: %d patients\n", length(high_risk_corrupted_idx)))
cat(sprintf("  * Identical allocations:     %d patients\n", intersection_ablation))
cat("====================================================\n\n")

# Step E: Generate Ablation Density Plot
ablation_df <- tibble(
  score = c(pred_corrupted, pred_clean),
  state = rep(c("Legacy State (Missing Month Data)", "Corrected State (Active Month Data)"), each = length(pred_clean))
)

p_ablation <- ggplot(ablation_df, aes(x = score, fill = state)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("#e74c3c", "#2ecc71")) +
  labs(
    title = "Date Feature Impact on Risk Stratification",
    subtitle = "Isolates the predictive shift of resolving the appt_month_num NA bug",
    x = "Calibrated DNA Risk Probability",
    y = "Density",
    fill = "Date Feature State"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave("outputs/audit/month_feature_ablation_density.png", plot = p_ablation, width = 8, height = 5)
message("Ablation audit completed successfully! Density plots saved to 'outputs/audit/'")
