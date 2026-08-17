# -------------------------------------------------------------------------
# 1. Setup & Configuration
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()

model_tune_results <- local({
  # Load prepped data and the recipe (decoupled from the prep script)
  train_data <- readRDS("data/processed/train_engineered.rds")
  dna_recipe <- readRDS("data/processed/dna_recipe.rds")
  
  # -------------------------------------------------------------------------
  # 2. Model Specification & Grid Setup
  # -------------------------------------------------------------------------
  rf_spec <- rand_forest(mtry = tune(),
                         trees = tune(),
                         min_n = tune()) %>%
    set_engine("ranger",
               num.threads = !!conf$num_threads,
               importance = "permutation") %>%
    set_mode("classification")
  
  # Juice the recipe once to find out how many columns we have for mtry()
  prepped_features <-
    prep(dna_recipe) %>% juice() %>% select(-all_of(conf$target_col))
  
  rf_grid <- grid_space_filling(
    mtry() %>% finalize(prepped_features),
    min_n(),
    trees(range = conf$tree_range),
    size = conf$grid_size
  )
  
  set.seed(123)
  dna_folds <-
    group_vfold_cv(train_data, v = conf$cv_folds, group = dim_patient_id)
  
  # -------------------------------------------------------------------------
  # 3. Parallel Execution & Tuning
  # -------------------------------------------------------------------------
  registerDoFuture()
  plan(
    multisession,
    workers = conf$num_workers,
    maxSizeOfObjects = 2000 * 1024 ^ 2
  )
  
  tic("Tidymodels grid tuning")
  
  fits <- workflow() %>%
    add_recipe(dna_recipe) %>%
    add_model(rf_spec) %>%
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