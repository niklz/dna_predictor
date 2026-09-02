# =========================================================================
# SCRIPT 04: FIT & CALIBRATE FINAL MODEL (CRASH-SAFE VERSION)
# =========================================================================
source("00_libraries_and_utils.R")
conf <- config::get()

train_data <- readRDS("data/processed/train_engineered.rds")

# deprecated
# if(conf$tune_model) {
#   source("02_tune_model.R")
# } else {
#   model_tune_results <- readRDS("data/processed/rf_tune.rds")
# }

# Fit the best model and calibrate on out-of-fold predictions
model <- local({
  
  params <- model_tune_results$best_params
  # -------------------------------------------------------------------------
  # 1. Reconstruct Untuned Workflow Blueprint [cite: 2]
  # -------------------------------------------------------------------------
  data_template <- head(train_data, 0)
  
  tuning_recipe <- build_trial_recipe(
    data_template = data_template,
    target_col    = conf$target_col,
    fct_other_prp = conf$fct_other_prp
  )
  
  rf_spec <- rand_forest(
    mtry  = tune(),
    trees = tune(),
    min_n = tune()
  ) %>%
    set_engine(
      "ranger",
      num.threads = 1,     # CRITICAL: Keeps background workers locked to 1 thread! [cite: 666]
      importance  = "none" # Keeps resample fitting incredibly lean [cite: 631]
    ) %>%
    set_mode("classification")
  
  wf <- workflow() %>%
    add_recipe(tuning_recipe) %>%
    add_model(rf_spec)
  
  # -------------------------------------------------------------------------
  # 2. Extract Unbiased Out-Of-Fold Predictions for Calibrator [cite: 415]
  # -------------------------------------------------------------------------
  message("Generating out-of-fold predictions for calibration...")
  
  # Recreate the folds for the best-model cross-validation run [cite: 3]
  set.seed(123)
  dna_folds <- rsample::group_vfold_cv(train_data, v = conf$cv_folds, group = "dim_patient_id")
  
  # Finalise blueprint with your winning parameters [cite: 406]
  best_wf <- wf %>% finalize_workflow(params)
  
  # ACTIVATE WORKSTATION PARALLEL CORES SAFELY [cite: 4]
  library(doFuture)
  registerDoFuture()
  plan(multisession, workers = conf$num_workers) # Launches parallel cores cleanly [cite: 4]
  
  # Run a highly parallel, thread-locked CV fit to harvest calibration predictions [cite: 413]
  cv_results <- fit_resamples(
    best_wf,
    resamples = dna_folds,
    metrics   = metric_set(pr_auc, roc_auc),
    control   = control_resamples(
      save_pred     = TRUE,       # Save predictions strictly for the winner [cite: 413]
      parallel_over = "resamples" # Parallelise strictly over folds [cite: 4]
    )
  )
  
  # Release parallel workers back to the system immediately [cite: 4]
  plan(sequential) 
  
  preds <- collect_predictions(cv_results)# Unbiased predictions harvested
  
  cal_model <- cal_estimate_logistic(
    preds,
    truth       = !!sym(conf$target_col),
    estimate    = c(.pred_DNA, .pred_attended),
    event_level = "first"
  )
  
  # -------------------------------------------------------------------------
  # 3. Fit Final Production Model (Multi-threading Override)
  # -------------------------------------------------------------------------
  message("Fitting final model on full cohort...")
  
  # Update parsnip engine to use maximum threads for the final fit 
  updated_spec <- rand_forest(
    mtry  = params$mtry,
    trees = params$trees,
    min_n = params$min_n
  ) %>%
    set_engine(
      "ranger",
      num.threads = parallel::detectCores() - 1, # Max threads for final fit
      importance  = "permutation"                # Safe to calculate now!
    ) %>%
    set_mode("classification")
  
  final_wf <- workflow() %>%
    add_recipe(tuning_recipe) %>%
    add_model(updated_spec)
  
  # Fit model on the full 1-million-row cohort sequentially [cite: 414]
  fit <- final_wf %>%
    fit(data = train_data)
  
  list(
    model      = fit,
    calibrator = cal_model,
    model_ver  = conf$model_ver
  )
})

saveRDS(model, "data/processed/rf_final_model.rds")
message("Success! Calibrated final model bundle saved.")