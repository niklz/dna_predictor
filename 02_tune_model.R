# -------------------------------------------------------------------------
# 1. Setup & Configuration
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()

model_tune_results <- local({
  # Load the engineered training data
  train_data <- readRDS("data/processed/train_engineered.rds")
  
  # -------------------------------------------------------------------------
  # 2. Define Workflow Builder (Prevents Environment Leakage)
  # -------------------------------------------------------------------------
  # Creating the workflow inside a temporary function isolates the environment.
  # Passing a 0-row template keeps the recipe object extremely lightweight.
  build_tuning_workflow <- function(target_col, data_template, fct_other_prp) {
    
    formula_obj <- as.formula(paste(target_col, "~ ."))
    # Strip any parent environment references from the formula
    environment(formula_obj) <- baseenv() 
    
    # Build a lightweight recipe using only the column blueprint
    dna_recipe <- recipe(formula_obj, data = data_template) %>%
      update_role(dim_patient_id, new_role = "id") %>%
      step_novel(all_nominal_predictors()) %>%
      step_unknown(all_nominal_predictors(), -imd) %>%
      step_other(all_nominal_predictors(), threshold = fct_other_prp) %>%
      step_nzv(all_predictors()) %>%
      step_impute_median(all_numeric_predictors())
    
    rf_spec <- rand_forest(
      mtry = tune(),
      trees = tune(),
      min_n = tune()
    ) %>%
      set_engine("ranger", importance = "permutation") %>%
      set_mode("classification")
    
    workflow() %>%
      add_recipe(dna_recipe) %>%
      add_model(rf_spec)
  }
  
  # Create a 0-row template containing only column names and types
  data_template <- head(train_data, 0)
  
  # Build the ultra-lightweight workflow container (typically < 100 KB)
  tuning_workflow <- build_tuning_workflow(
    target_col    = conf$target_col,
    data_template = data_template,
    fct_other_prp = conf$fct_other_prp
  )
  
  # -------------------------------------------------------------------------
  # 3. Create Resamples and Parallel Grid Tuning
  # -------------------------------------------------------------------------
  set.seed(123)
  dna_folds <- group_vfold_cv(train_data, v = conf$cv_folds, group = dim_patient_id)
  
  # Set up the search grid (juice the recipe locally once)
  prepped_features <- prep(readRDS("data/processed/dna_recipe.rds")) %>% 
    juice() %>% 
    select(-all_of(conf$target_col))
  
  rf_grid <- grid_space_filling(
    mtry() %>% finalize(prepped_features),
    min_n(),
    trees(range = conf$tree_range),
    size = conf$grid_size
  )
  
  # Set up parallel execution with safety parameters
  registerDoFuture()
  plan(
    multisession, 
    workers = conf$num_workers,
    # This option is now safe because the workflow size has dropped dramatically
    future.globals.maxSize = 1000 * 1024^2 
  )
  
  tic("Tidymodels grid tuning")
  
  fits <- tuning_workflow %>%
    tune_grid(
      resamples = dna_folds,
      grid = rf_grid,
      metrics = metric_set(pr_auc, roc_auc),
      control = control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        parallel_over = "everything"
      )
    )
  
  toc()
  
  # -------------------------------------------------------------------------
  # 4. Diagnostics & Plotting
  # -------------------------------------------------------------------------
  best_params <- fits %>% select_best(metric = "pr_auc")
  
  hyper_param_trace <- fits %>%
    collect_metrics() %>%
    filter(.metric == "pr_auc") %>%
    select(mean, mtry:min_n) %>%
    pivot_longer(mtry:min_n, values_to = "value", names_to = "parameter") %>%
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(alpha = 0.8, show.legend = FALSE) +
    facet_wrap(~ parameter, scales = "free_x") +
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
  
  # Performance plots using the best params
  roc_plot <-
    collect_predictions(fits, parameters = best_params) |> (\(preds) {
      auc_val <-
        preds |> roc_auc(truth = !!sym(conf$target_col), .pred_DNA) |> pull(.estimate)
      preds |>
        roc_curve(truth = !!sym(conf$target_col), .pred_DNA) |>
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_abline(slope = 1,
                    linetype = 2,
                    alpha = 0.4) +
        geom_path(linewidth = 1, color = "midnightblue") +
        coord_equal() +
        theme_bw() +
        labs(title = "ROC curve", subtitle = paste0("ROC AUC: ", round(auc_val, 3)))
    })()
  
  pr_plot <-
    collect_predictions(fits, parameters = best_params) |> (\(preds) {
      pr_val <-
        preds |> pr_auc(truth = !!sym(conf$target_col), .pred_DNA) |> pull(.estimate)
      preds |>
        pr_curve(truth = !!sym(conf$target_col), .pred_DNA) %>%
        ggplot(aes(x = recall, y = precision)) +
        geom_path(linewidth = 1, color = "midnightblue") +
        geom_hline(yintercept = 0.05,
                   lty = 2,
                   color = "red") +
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
    roc_plot = roc_plot
  )
  
  model_tune_results
})


saveRDS(model_tune_results, "data/processed/rf_tune.rds")

message("Tuning complete. Diagnostics saved to data/processed/")