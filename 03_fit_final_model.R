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
  
  model_tune_results <- readRDS("data/processed/rf_tune.rds")
  params <- model_tune_results$best_params
  # -------------------------------------------------------------------------
  # 1. Reconstruct untuned workflow
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
      num.threads = 1,
      importance  = "none" # Keeps resample fitting lean
    ) %>%
    set_mode("classification")
  
  wf <- workflow() %>%
    add_recipe(tuning_recipe) %>%
    add_model(rf_spec)
  
  # -------------------------------------------------------------------------
  # 2. Extract unbiased out-of-fold predictions for calibrator
  # -------------------------------------------------------------------------
  message("Generating out-of-fold predictions for calibration...")
  
  # Recreate the folds for the best-model cross-validation run
  set.seed(123)
  dna_folds <- rsample::group_vfold_cv(train_data, v = conf$cv_folds, group = "dim_patient_id")
  
  # Finalise blueprint with winning parameters 
  best_wf <- wf %>% finalize_workflow(params)
  
  # LAPTOP-SAFE SEQUENTIAL RUN
  library(doFuture)
  registerDoFuture()
  plan(sequential)
  
  # Run the CV fit to harvest calibration predictions
  cv_results <- fit_resamples(
    best_wf,
    resamples = dna_folds,
    metrics   = metric_set(pr_auc, roc_auc),
    control   = control_resamples(
      save_pred     = TRUE,
      parallel_over = "resamples"
    )
  )
  

  preds <- collect_predictions(cv_results)
  
  cal_model <- cal_estimate_logistic(
    preds,
    truth       = !!sym(conf$target_col),
    estimate    = c(.pred_DNA, .pred_attended),
    event_level = "first"
  )
  
  # -------------------------------------------------------------------------
  # 3. Fit final production model
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
      importance  = "permutation"                # Importance only computed for final fit
    ) %>%
    set_mode("classification")
  
  final_wf <- workflow() %>%
    add_recipe(tuning_recipe) %>%
    add_model(updated_spec)
    
  fit <- final_wf %>%
    fit(data = train_data)
  
  list(
    model      = fit,
    calibrator = cal_model,
    predictions = preds,
    model_ver  = conf$model_ver
  )
})

saveRDS(model, "data/processed/rf_final_model.rds")
message("Success! Calibrated final model bundle saved.")