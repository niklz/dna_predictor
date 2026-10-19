# =========================================================================
# SCRIPT 07: LEARNING CURVE & MODEL SIZE OPTIMISER
# =========================================================================
# This script programmatically evaluates model performance (PR AUC, ROC AUC),
# training time, and physical file size across varying training set sizes
# (e.g., 10k, 50k, 100k, 250k, 500k rows) to find the "sweet spot".

library(dplyr)
library(tidymodels)
library(tictoc)
library(readr)

# Load configuration and utilities
source("00_libraries_and_utils.R")
conf <- config::get()

# Set up logging directories
dir.create("outputs/diagnostics", showWarnings = FALSE, recursive = TRUE)

# 1. Load the raw engineered training data (assumed to be ~1M rows)
message("Loading full training cohort...")
train_full <- readRDS("data/processed/train_engineered.rds")
total_rows <- nrow(train_full)
message(sprintf("Loaded %s rows of training data.", format(total_rows, big.mark=",")))

# Define sub-sampling sizes to test
sampling_sizes <- c(10000, 50000, 100000, 250000, 500000)
if (total_rows > 500000) {
  sampling_sizes <- c(sampling_sizes, total_rows)
}

# Pre-allocate results dataframe
curve_results <- data.frame(
  sample_size       = numeric(),
  pct_of_data       = numeric(),
  training_time_sec = numeric(),
  file_size_mb      = numeric(),
  pr_auc_val        = numeric(),
  roc_auc_val       = numeric(),
  stringsAsFactors  = FALSE
)

# 2. Extract best parameters
model_tune  = readRDS("data/processed/rf_tune.rds")
best_params = model_tune$best_params

# Define model specification (ensuring num.threads = 1 to prevent thrashing on loops)
rf_spec <- rand_forest(
  mtry  = best_params$mtry,
  trees = best_params$trees,
  min_n = best_params$min_n
) %>%
  set_engine("ranger", num.threads = conf$num_threads, importance = "none") %>%
  set_mode("classification")

# Define recipe template on empty schema
data_template <- head(train_full, 0)
base_recipe   <- build_trial_recipe(data_template, target_col = conf$target_col, fct_other_prp = conf$fct_other_prp)

# Create a validation/test set to evaluate out-of-sample performance
set.seed(42)
split_idx <- rsample::initial_split(train_full, prop = 0.8, strata = !!sym(conf$target_col))
train_pool <- rsample::training(split_idx)
eval_set   <- rsample::testing(split_idx)

# Pre-prep evaluation features
eval_engineered <- eval_set

# 3. Loop through sample sizes
for (n in sampling_sizes) {
  message(sprintf("\n--- Evaluating Sample Size: %s rows ---", format(n, big.mark=",")))
  
  # Stratified sub-sampling to preserve class balance
  set.seed(123)
  sample_data <- train_pool %>%
    group_by(!!sym(conf$target_col)) %>%
    slice_sample(prop = pmin(1, n / nrow(train_pool)), replace = FALSE) %>%
    ungroup()
  
  # Time the model training
  tic()
  
  # Build workflow and fit
  wf <- workflow() %>%
    add_recipe(base_recipe) %>%
    add_model(rf_spec)
  
  fitted_model <- wf %>% fit(data = sample_data)
  
  time_elapsed <- toc(quiet = TRUE)
  train_time <- time_elapsed$toc - time_elapsed$tic
  
  # Test serialization file size (saving temporarily without compression to check raw size)
  temp_file <- tempfile(fileext = ".rds")
  # Strip environments first to make a fair comparison
  library(butcher)
  cleaned_fit <- butcher::butcher(fitted_model)
  
  saveRDS(cleaned_fit, temp_file, compress = FALSE)
  file_sz <- file.size(temp_file) / (1024^2) # Size in MB
  unlink(temp_file)
  
  # Generate out-of-sample predictions
  eval_preds <- predict(cleaned_fit, new_data = eval_engineered, type = "prob") %>%
    bind_cols(eval_engineered)
  
  # Compute performance metrics
  pr_res <- eval_preds %>% pr_auc(truth = !!sym(conf$target_col), .pred_DNA)
  roc_res <- eval_preds %>% roc_auc(truth = !!sym(conf$target_col), .pred_DNA)
  
  pr_auc_val  <- pr_res$.estimate
  roc_auc_val <- roc_res$.estimate
  
  # Append to results
  curve_results <- rbind(curve_results, data.frame(
    sample_size       = n,
    pct_of_data       = (n / total_rows) * 100,
    training_time_sec = round(train_time, 2),
    file_size_mb      = round(file_sz, 2),
    pr_auc_val        = round(pr_auc_val, 4),
    roc_auc_val       = round(roc_auc_val, 4)
  ))
  
  # Print iteration metrics
  cat(sprintf("Training Time: %.2f sec | File Size: %.1f MB | PR AUC: %.4f | ROC AUC: %.4f\n", 
              train_time, file_sz, pr_auc_val, roc_auc_val))
}

# 4. Save results report
write_csv(curve_results, "outputs/diagnostics/learning_curve_summary.csv")
message("\nSuccess! Learning curve completed. Summary saved to outputs/diagnostics/learning_curve_summary.csv")

# 5. Output programmatic advice on the optimal size
sweet_spot <- curve_results %>%
  filter(pr_auc_val >= (max(pr_auc_val) - 0.002)) %>% # Accept up to 0.2% loss from absolute max
  arrange(sample_size) %>%
  slice(1)

cat("\n====================================================\n")
cat("      LEARNING CURVE SWEET-SPOT ANALYSIS            \n")
cat("====================================================\n")
cat(sprintf("Absolute Best Performance: PR AUC = %.4f (N = %s)\n", 
            max(curve_results$pr_auc_val), format(max(curve_results$sample_size), big.mark=",")))
cat(sprintf("Recommended Sweet Spot:   N = %s rows (%.1f%% of dataset)\n", 
            format(sweet_spot$sample_size, big.mark=","), sweet_spot$pct_of_data))
cat(sprintf("-> Saves %.1f%% in training time (%.1f vs %.1f seconds)\n", 
            (1 - sweet_spot$training_time_sec / max(curve_results$training_time_sec)) * 100,
            sweet_spot$training_time_sec, max(curve_results$training_time_sec)))
cat(sprintf("-> Shrinks model file size by %.1f%% (%.1f vs %.1f MB)\n", 
            (1 - sweet_spot$file_size_mb / max(curve_results$file_size_mb)) * 100,
            sweet_spot$file_size_mb, max(curve_results$file_size_mb)))
cat("====================================================\n\n")
