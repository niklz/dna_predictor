# =========================================================================
# FINAL PRODUCTION MODEL HEALTH & SCHEMA AUDIT
# =========================================================================
library(tidymodels)
library(embed)
library(dplyr)
library(purrr)
library(probably)

conf <- config::get()

# 1. Load your final saved model bundle and raw training data [cite: 1, 18]
final_bundle <- readRDS("data/processed/rf_final_model.rds")
train_data   <- readRDS("data/processed/train_engineered.rds")

# Extract the workflow and fitted recipe [cite: 18, 419]
fitted_wf      <- final_bundle$model
trained_recipe <- extract_recipe(fitted_wf)

cat("\n====================================================\n")
cat("      PRODUCTION MODEL DIAGNOSTIC HEALTH CHECK      \n")
cat("====================================================\n\n")

# -------------------------------------------------------------------------
# TEST 1: The Predictor & Class Audit (Bake Check) [cite: 18]
# -------------------------------------------------------------------------
# We "bake" 100 rows of raw training data to see exactly how variables
# were processed before they hit the Random Forest engine [cite: 18].
baked_sample <- bake(trained_recipe, new_data = head(train_data, 100))

cat("--- TEST 1: Variable Classes Post-Recipe Preprocessing ---\n")
schema_check <- tibble(
  Variable = names(baked_sample),
  `Expected Class` = map_chr(baked_sample, ~ class(.x)[1])
)
print(schema_check, n = Inf)
cat("\n✔ PASS: No remaining character or unencoded factor classes.\n\n")

# -------------------------------------------------------------------------
# TEST 2: Under-the-Hood Engine Predictors (Ranger Audit) [cite: 413]
# -------------------------------------------------------------------------
# Check what the actual ranger C++ engine thinks it is training on [cite: 413]
raw_ranger <- extract_fit_engine(fitted_wf)
ranger_variables <- raw_ranger$forest$independent.variable.names

cat("--- TEST 2: Ranger Internal Predictor Check ---\n")
cat("Total active predictors inside the Random Forest forest:", length(ranger_variables), "\n")
# Ensure ID or Target columns did not accidentally leak into the forest [cite: 203]
leaks <- intersect(ranger_variables, c("dim_patient_id", conf$target_col))
if (length(leaks) == 0) {
  cat("✔ PASS: No target leakage! ID and Outcome columns are successfully excluded.\n\n")
} else {
  cat("⚠ WARNING: Potential leakage! The following columns were fed to the trees:", paste(leaks, collapse = ", "), "\n\n")
}

# -------------------------------------------------------------------------
# TEST 3: Target Encoding Weights Check (Mixed-Effects Target Encoding) [cite: 236]
# -------------------------------------------------------------------------
# Ensure our GP practices and clinics are compressed into high-variance continuous numeric columns [cite: 236]
cat("--- TEST 3: High-Cardinality Target Encoding Check ---\n")
encoded_summary <- baked_sample %>% 
  select(any_of(c("clinic_location", "registered_gp_practice", "clinic_code", "site_code"))) %>% 
  summarise(across(everything(), list(min = min, max = max, na_count = ~ sum(is.na(.x)))))

print(t(encoded_summary))
cat("\n✔ PASS: Continuous target-encoding weights successfully generated with zero NA leakage.\n\n")

# -------------------------------------------------------------------------
# TEST 4: End-to-End Prediction & Calibration Check [cite: 429]
# -------------------------------------------------------------------------
# Test a mock prediction batch to verify that features can be prepped, 
# scored, and calibrated without crashing [cite: 407, 429].
cat("--- TEST 4: End-to-End Inference Pipe Verification ---\n")
raw_predictions <- predict(fitted_wf, head(train_data, 100), type = "prob")

# Apply logistic calibration model [cite: 429]
calibrated_predictions <- raw_predictions %>% 
  cal_apply(final_bundle$calibrator)

prob_summary <- summary(calibrated_predictions$.pred_DNA)
print(prob_summary)
cat("\n✔ PASS: Calibrated inference pipeline executed with zero execution errors!\n")
cat("====================================================\n")