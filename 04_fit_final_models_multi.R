# =========================================================================
# SCRIPT 04: FIT & CALIBRATE FINAL MODELS (MULTI-VARIANT DESKTOP VERSION)
# =========================================================================
# Sourced for desktop execution to leverage Core i9-9980XE (36 logical threads).
# Fits, calibrates, and outputs dual-package RDS bundles (Full and Lean)
# across multiple training cohort sizes: 100k, 250k, 500k, and Full (944k).

source("00_libraries_and_utils.R")
conf <- config::get()

# Set up logging directories
dir.create("data/processed/models", showWarnings = FALSE, recursive = TRUE)

# 1. Load the full training data
message("Loading full training cohort...")
train_full <- readRDS("data/processed/train_engineered.rds")
total_rows <- nrow(train_full)
message(sprintf("Loaded %s rows of training data.", format(total_rows, big.mark=",")))

# 2. Extract best parameters from existing tuning grid
message("Loading tuning hyperparameters...")
model_tune_results <- readRDS("data/processed/rf_tune.rds")
best_params        <- model_tune_results$best_params
message(sprintf("Using optimal tuned parameters: mtry = %d, trees = %d, min_n = %d", 
                best_params$mtry, best_params$trees, best_params$min_n))

# -------------------------------------------------------------------------
# Dynamic Scoping Environment Stripper (To prevent lexical closures)
# -------------------------------------------------------------------------
strip_lexical_envs <- function(x) {
  if (is.list(x) && !is.object(x)) {
    x <- lapply(x, strip_lexical_envs)
  }
  if (inherits(x, "formula") || inherits(x, "terms")) {
    environment(x) <- baseenv()
  }
  if (!is.null(attr(x, ".Environment"))) {
    attr(x, ".Environment") <- baseenv()
  }
  return(x)
}

# Define training cohort sizes to build
cohort_sizes <- c(500000, total_rows)

# -------------------------------------------------------------------------
# LOOP THROUGH COHORT SIZES
# -------------------------------------------------------------------------
for (n in cohort_sizes) {
  size_label <- if (n == total_rows) "full" else sprintf("%dk", n / 1000)
  message(sprintf("\n===================================================="))
  message(sprintf("   TRAINING COHORT VARIANT: %s rows (%s)", format(n, big.mark=","), size_label))
  message(sprintf("===================================================="))
  
  # A. Stratified sub-sampling to preserve class balance (DNA rate)
  set.seed(123)
  train_subset <- if (n == total_rows) {
    train_full
  } else {
    train_full %>%
      group_by(!!sym(conf$target_col)) %>%
      slice_sample(prop = pmin(1, n / total_rows), replace = FALSE) %>%
      ungroup()
  }
  
  # B. Reconstruct Workflow Blueprints
  data_template <- head(train_subset, 0)
  tuning_recipe <- build_trial_recipe(
    data_template = data_template, 
    target_col    = conf$target_col,
    fct_other_prp = conf$fct_other_prp
  )
  
  # Resample specification (Lock workers to 1 thread during CV to avoid thrashing)
  # FIXED: Keep parameters untuned (tune()) and finalize_workflow() below.
  # This completely bakes the actual parameters as constants into the workflow structure, 
  # ensuring parallel workers on Windows socket nodes never look for 'best_params' in global parent scopes.
  rf_cv_spec <- rand_forest(
    mtry  = tune(), 
    trees = tune(), 
    min_n = tune()
  ) %>% 
    set_engine("ranger", num.threads = 1, importance = "none") %>% 
    set_mode("classification")
  
  cv_wf <- workflow() %>% 
    add_recipe(tuning_recipe) %>% 
    add_model(rf_cv_spec) %>%
    finalize_workflow(best_params)
  
  # C. Generate 10-Fold CV Resamples
  # Note: Stratification is omitted in group cross-validation because 
  # dna_outcome is not constant within each patient ID (dim_patient_id) group.
  set.seed(123)
  dna_folds <- rsample::group_vfold_cv(
    train_subset, 
    v = conf$cv_folds, 
    group = "dim_patient_id"
  )
  
  # D. Parallel Evaluation of Folds
  message("Running parallel out-of-fold evaluations for calibration...")
  library(doFuture)
  registerDoFuture()
  # Desktop optimization: 10 parallel workers to match your 10 folds!
  plan(multisession, workers = 10, maxSizeOfObjects = 3000 * 1024^2)
  
  cv_results <- fit_resamples(
    cv_wf,
    resamples = dna_folds,
    metrics   = metric_set(pr_auc, roc_auc),
    control   = control_resamples(
      save_pred     = TRUE,
      parallel_over = "resamples"
    )
  )
  
  # Close parallel session to release workers
  plan(sequential)
  
  # E. Collect Unbiased Predictions and Build Calibrator
  preds <- collect_predictions(cv_results)
  cal_model <- cal_estimate_logistic(
    preds,
    truth       = !!sym(conf$target_col),
    estimate    = c(.pred_DNA, .pred_attended),
    event_level = "first"
  )
  
  # F. Fit Final Production Model (Multi-threading enabled)
  message("Fitting final model on complete subset...")
  final_spec <- rand_forest(
    mtry  = tune(),
    trees = tune(),
    min_n = tune()
  ) %>% 
    set_engine("ranger", num.threads = 16, importance = "permutation") %>% 
    set_mode("classification")
  
  final_wf <- workflow() %>% 
    add_recipe(tuning_recipe) %>% 
    add_model(final_spec) %>%
    finalize_workflow(best_params)
  
  fit <- final_wf %>% fit(data = train_subset)
  
  # G. Construct and Save Bundles
  message("Constructing bundles...")
  
  # --- PACKAGE 1: THE FULL EVALUATION BUNDLE ---
  # Retains predictions table (needed for threshold evaluation)
  # Stripped of lexical environment closures to prevent serialization traps
  full_bundle <- list(
    model       = fit,
    calibrator  = cal_model,
    predictions = tibble::as_tibble(as.data.frame(preds)), # Strip tidy evaluation metadata
    model_ver   = sprintf("%s_%s_full", conf$model_ver, size_label)
  )
  
  full_bundle <- strip_lexical_envs(full_bundle)
  full_bundle$model <- butcher::butcher(full_bundle$model)
  full_bundle$calibrator <- butcher::butcher(full_bundle$calibrator)
  
  full_path <- sprintf("data/processed/models/rf_final_model_%s_full.rds", size_label)
  saveRDS(full_bundle, full_path, compress = "gzip")
  
  # --- PACKAGE 2: THE LEAN PREDICTION BUNDLE ---
  # Stripped of predictions table (saves vast amounts of space)
  # Double-butchered and stripped of all environment bindings
  lean_bundle <- list(
    model       = full_bundle$model,
    calibrator  = full_bundle$calibrator,
    predictions = NULL, # Completely stripped!
    model_ver   = sprintf("%s_%s_lean", conf$model_ver, size_label)
  )
  
  # Save with zero compression for instant loadRDS() times on laptop serving pipelines
  lean_path <- sprintf("data/processed/models/rf_final_model_%s_lean.rds", size_label)
  saveRDS(lean_bundle, lean_path, compress = FALSE)
  
  # Print diagnostic sizes
  full_size <- file.size(full_path) / 1024^2
  lean_size <- file.size(lean_path) / 1024^2
  message(sprintf("-> [COMPLETE] Saved Full Bundle (%s): %.2f MB", size_label, full_size))
  message(sprintf("-> [COMPLETE] Saved Lean Bundle (%s): %.2f MB", size_label, lean_size))
}

message("\n====================================================")
message("   ALL MODELS TRAINED AND SERIALISED SUCCESSFULLY!")
message("====================================================")
