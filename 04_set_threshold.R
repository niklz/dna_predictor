# =========================================================================
# SCRIPT 05: SET OPTIMAL ML THRESHOLDS
# =========================================================================
# Uses cross-validated predictions to find a threshold that satisfies both
# statistical trial requirements and operational coordinator capacity.

source("00_libraries_and_utils.R")

optimal_bounds <- local({

  # 1. Load tuning results (which contains the hold-out predictions)
  model_tune_results <- readRDS("data/processed/rf_tune.rds")
  model_final <- readRDS("data/processed/rf_final_model.rds")
  best_params <- model_tune_results$best_params
  
  # 2. Extract out-of-fold predictions from the best CV model & apply calibration
  # This provides an unbiased estimate of the historical score distribution
  cv_preds <- model_final$predictions
  
  # 3. Generate PR Curve data and raw historical scores
  pr_curve_data <-
    cv_preds %>% pr_curve(truth = dna_outcome, .pred_DNA)
  historical_scores <- cv_preds %>% pull(.pred_DNA)
  
  # 4. Calculate Optimal Thresholds
  # Using your clinic's specific operational parameters
  optimal_bounds <- find_optimal_ml_thresholds(
    pr_curve_data     = pr_curve_data,
    historical_scores = historical_scores,
    v_week            = 308,
    # Total weekly clinic volume
    w                 = 26,
    # Trial duration in weeks
    n_req             = 800,
    # Minimum trial sample size
    a_int             = 0.50,
    # 50/50 Intervention split
    r_call            = 0.70,
    # 70% need a manual call (no response to auto text)
    c_max             = 30    # Coordinator call limit per week
  )
  
  # 5. Lock in the Final Production Threshold
  # If feasible, we default to the ceiling threshold (most conservative for staff
  # while still hitting power). If invalid, we enforce the ceiling anyway so
  # operations don't crash, and flag a warning.
  final_production_threshold <- optimal_bounds$ceiling_threshold
  
  if (!optimal_bounds$feasible) {
    warning(
      "Safe zone invalid! Enforcing ceiling threshold to protect coordinator capacity. Trial may run longer than `w` weeks to reach `n_req`."
    )
  }
  optimal_bounds
})



# 6. Save the threshold artifact for the Inference script
saveRDS(
  list(
    bounds = optimal_bounds,
    production_threshold = optimal_bounds$ceiling_threshold
  ),
  file = "data/processed/risk_threshold.RDS"
)

print(paste(
  "Final risk threshold locked at >=",
  round(optimal_bounds$ceiling_threshold, 4)
))
