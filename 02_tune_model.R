# -------------------------------------------------------------------------
# 1. Setup & Configuration
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()


# # Single sequential fold timing test
# 
# local({
#   train_data <- readRDS("data/processed/train_engineered_tune.rds")
#   
#   # 2. Lock SPEC to single thread
#   rf_spec <- rand_forest(mtry = 3, trees = 500, min_n = 6) %>%
#     set_engine("ranger", num.threads = 1, importance = "none") %>%
#     set_mode("classification")
#   
#   # 3. Create a 1-fold rset using safe row-bracket indexing!
#   set.seed(123)
#   all_folds <- group_vfold_cv(train_data, v = 10, group = "dim_patient_id")
#   
#   single_fold <- rsample::manual_rset(all_folds$splits[1], ids = "Fold1")
#   
#   # 4. Bind and run sequentially (no parallel setup!)
#   wf <- workflow() %>% 
#     add_recipe(readRDS("data/processed/dna_recipe.rds")) %>% 
#     add_model(rf_spec)
#   
#   library(tictoc)
#   tic("Single sequential fold")
#   fit_resamples(
#     wf,
#     resamples = single_fold,
#     metrics   = metric_set(pr_auc, roc_auc),
#     control   = control_resamples(save_pred = FALSE)
#   )
#   toc()
# })



model_tune_results <- local({
  # Load the engineered training data
  train_data <- readRDS("data/processed/train_engineered_tune.rds") 
  
  rf_spec <- rand_forest(mtry = tune(),
                         trees = tune(),
                         min_n = tune()) %>%
    set_engine(
      "ranger",
      num.threads = 1,          
      importance  = "none"      
      ) %>%
    set_mode("classification")
  
  
  # Create a 0-row template containing only column names and types
  data_template <- head(train_data, 0)
  
  # Build the recipe using the exact same centralized configuration [cite: 11]
  tuning_recipe <- build_trial_recipe(
    data_template = data_template,
    target_col    = conf$target_col,
    fct_other_prp = conf$fct_other_prp
  )
  
  # Bind the recipe and model specification to your workflow
  tuning_workflow <- workflow() %>%
    add_recipe(tuning_recipe) %>%
    add_model(rf_spec)
  
  # -------------------------------------------------------------------------
  # 3. Create Resamples and Parallel Grid Tuning
  # -------------------------------------------------------------------------
  set.seed(123)
  dna_folds <- group_vfold_cv(
    train_data,
    v = conf$cv_folds,
    group = "dim_patient_id"
    )

  
  # Set up the search grid (juice the recipe locally once)
  prepped_features <- prep(readRDS("data/processed/dna_recipe.rds"), training = train_data) %>%
    juice() %>%
    select(-all_of(conf$target_col))
  
  rf_grid <- grid_space_filling(
    mtry() %>% finalize(prepped_features),
    min_n(),
    trees(range = conf$tree_range),
    size = conf$grid_size
  )
  
  # Set up parallel execution with safety parameters
  options(future.globals.maxSize = 8000 * 1024^2)
  registerDoFuture()
  plan(multisession, workers = conf$num_workers)
  
  tic("Tidymodels grid tuning")
  
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  control_safe <- control_grid(
    save_pred     = FALSE,       # Discard intermediate predictions
    save_workflow = FALSE,       # Discard bloated workflows to prevent memory leaks
    parallel_over = "resamples", # Parallelise strictly over folds (reduces socket copies from 250 to 10)
    verbose       = TRUE
  )
  
    
    fits <- tuning_workflow %>%
      tune_grid(
        resamples = dna_folds,
        grid = rf_grid,
        metrics = metric_set(pr_auc, roc_auc),
        control = control_safe
      )
  
  toc()
  
  # -------------------------------------------------------------------------
  # 4. Diagnostics & Plotting (Crash-Safe Edition)
  # -------------------------------------------------------------------------
  best_params <- fits %>% select_best(metric = "pr_auc")
  
  # A. Hyperparameter trace and parallel coordinates (These are safe with collect_metrics!) [cite: 4, 5]
  hyper_param_trace <- fits %>%
    collect_metrics() %>%
    filter(.metric == "pr_auc") %>%
    select(mean, mtry:min_n) %>%
    pivot_longer(mtry:min_n, values_to = "value", names_to = "parameter") %>%
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(alpha = 0.8, show.legend = FALSE) +
    facet_wrap( ~ parameter, scales = "free_x") +
    labs(x = NULL, y = "PR AUC")
  
  para_coord_plot <- fits %>%
    collect_metrics() %>%
    filter(.metric == "pr_auc") %>%
    select(mtry, trees, min_n, mean) %>%
    ggparcoord(
      columns = 1:3,
      groupColumn = 4,
      scale = "uniminmax",
      showPoints = TRUE,
      alphaLines = 0.7
    ) +
    scale_color_viridis_c(option = "viridis", name = "Mean PR-AUC") +
    theme_minimal() +
    labs(
      title = "Parallel Coordinates Plot of Random Forest Tuning",
      subtitle = "Higher PR-AUC paths are highlighted in yellow/green",
      x = "Hyperparameters",
      y = "Normalized Scale (0 to 1)"
    )
  
  # -------------------------------------------------------------------------
  # B. RE-FIT THE WINNER TO HARVEST UNBIASED PREDICTIONS [cite: 413]
  # -------------------------------------------------------------------------
  # This runs in seconds and is 100% memory-safe!
  message("Generating out-of-fold predictions for the best model...")
  
  best_workflow <- tuning_workflow %>% 
    finalize_workflow(best_params)
  
  best_cv_results <- fit_resamples(
    best_workflow,
    resamples = dna_folds,
    metrics   = metric_set(pr_auc, roc_auc),
    control   = control_resamples(save_pred = TRUE) # Only save predictions for the winner! [cite: 413]
  )
  
  best_preds <- collect_predictions(best_cv_results) # Unbiased predictions harvested! [cite: 413, 415]
  
  # -------------------------------------------------------------------------
  # C. Generate ROC & PR Plots (Using the safe best_preds) [cite: 6]
  # -------------------------------------------------------------------------
  roc_plot <- best_preds |> (\(preds) {
    auc_val <- preds |> roc_auc(truth = !!sym(conf$target_col), .pred_DNA) |> pull(.estimate)
    preds |>
      roc_curve(truth = !!sym(conf$target_col), .pred_DNA) |>
      ggplot(aes(x = 1 - specificity, y = sensitivity)) +
      geom_abline(slope = 1, linetype = 2, alpha = 0.4) +
      geom_path(linewidth = 1, color = "midnightblue") +
      coord_equal() +
      theme_bw() +
      labs(title = "ROC curve", subtitle = paste0("ROC AUC: ", round(auc_val, 3)))
  })()
  
  pr_plot <- best_preds |> (\(preds) {
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
  
  
  # -------------------------------------------------------------------------
  # 5. Save Tuning Output
  # -------------------------------------------------------------------------
  model_tune_results <- list(
    tune_res = fits,
    best_params = best_params,
    param_trace = hyper_param_trace,
    para_coord_plot = para_coord_plot,
    pr_plot = pr_plot,
    roc_plot = roc_plot,
    model_ver = conf$model_ver
  )
  
  model_tune_results
})


saveRDS(model_tune_results, "data/processed/rf_tune.rds")
saveRDS(
  model_tune_results,
  "S:/Finance/Shared Area/BNSSG - BI/8 Modelling and Analytics/working/nh/projects/dna_predictor/data/rf_tune.RDS"
)

message("Tuning complete. Diagnostics saved to data/processed/")