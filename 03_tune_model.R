source("02_data_prep.R")


model_tune <- local({
  
  collect_ethnicity <- function(x) {
    dplyr::case_when(
      x %in% c("Not Stated", "Not Known", "Not Collected At This Time", "Not Set") ~ "Unknown",
      x %in% c("White - British", "White - Any other White background", "White - Irish", "Any other White background") ~ "White",
      x %in% c("Asian or Asian British - Indian", "Asian - Indian or British Indian", 
               "Asian or Asian British - Pakistani", "Asian - Pakistani or British Pakistani",
               "Asian or Asian British - Bangladeshi", "Asian - Bangladeshi or British Bangladeshi",
               "Asian or Asian British - Any other Asian background", "Asian - Any other Asian background") ~ "Asian",
      x %in% c("Black or Black British - Caribbean", "Black - Caribbean or Black British Caribbean",
               "Black Caribbean", "Black or Black British - African", 
               "Black  - African or British African", "Black or Black British - Any other Black background",
               "Black - Any other Black background") ~ "Black",
      x %in% c("Mixed - White and Black African", "Mixed - Any other background", 
               "Mixed - White and Black Caribbean", "Mixed - White and Asian", "Any other Mixed Background") ~ "Mixed",
      x %in% c("Other Ethnic Group - Chinese", "Chinese", "Any Other Ethnic Group", "Any other Ethnic Group") ~ "Other",
      TRUE ~ "Unknown"
    )
  }
  
  aggregate_ethnicity_high_level <- function(x) {
    # Standardise string inputs
    x_clean <- trimws(tolower(as.character(x)))
    
    dplyr::case_when(
      # Unknown / Missing / Unspecified
      x_clean %in% c("unknown", "not stated", "not known", 
                     "not collected at this time", "not set") ~ "Unknown",
      
      # White
      grepl("^white", x_clean) ~ "White",
      
      # Asian / Asian British
      grepl("^asian", x_clean) ~ "Asian",
      
      # Black / Black British
      grepl("^black", x_clean) ~ "Black",
      
      # Mixed / Multiple Ethnic Groups
      grepl("^mixed", x_clean) ~ "Mixed",
      
      # Chinese and Other Ethnic Groups
      grepl("chinese", x_clean) | grepl("other", x_clean) ~ "Other",
      
      # Catch-all fallback
      TRUE ~ "Unknown"
    )
  }
  
# --- 2. The Recipe (Pre-processing Pipeline) ---
# Tidymodels handles "knowledge separation" automatically.
  dna_recipe <- recipe(dna_outcome ~ ., data = train_raw) %>%
    update_role(dim_patient_id, new_role = "id") %>%
    step_mutate(
      ethnicity_clean = clean_ethnicity(ethnicity),
      ethnicity_group = aggregate_ethnicity_high_level(ethnicity_clean),
      # Set the baseline reference level to the most frequent category
      ethnicity_group = factor(ethnicity_group) %>% relevel(ref = "White"),
      appt_date = as.Date(substring(appt_month, 1, 10), format = "%d/%m/%Y"),
      appt_dow = factor(weekdays(appt_date)),
      appt_month_num = as.factor(format(appt_date, "%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      lead_time_days_log = log1p(pmax(0, lead_time_days)),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      appt_hour_sin = sin(2 * pi * appt_hour / 24),
      appt_hour_cos = cos(2 * pi * appt_hour / 24),
      has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
    ) %>%
    step_rm(appt_hour,
            lead_time_days,
            appt_date,
            appt_month,
            prev_dna_ly,
            ethnicity,
            ethnicity_clean) %>%
    step_novel(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors(), -imd) %>%
    step_other(all_nominal_predictors(), threshold = fct_other_prp) %>%
    step_zv(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_lencode_mixed(any_of(
      c(
        "clinic_location",
        "clinic_code",
        "site_code",
        "registered_gp_practice"
      )
    ), outcome = vars(dna_outcome)) %>%
    step_impute_median(all_numeric_predictors())
  
  # --- 3. Model Specification ---
  rf_spec <- rand_forest(mtry = tune(),
                         trees = tune(),
                         min_n = tune()) %>%
    set_engine(
      "ranger",
      num.threads = 4,
      #num.threads = parallel::detectCores(),
      importance = "permutation"
    ) %>%
    set_mode("classification")
  
  prepped_features <- prep(dna_recipe) %>% juice() %>% select(-dna_outcome)
  
  rf_grid <- grid_space_filling(mtry() %>% finalize(prepped_features),
                                min_n(),
                                trees(range = tree_range),
                                size = 25)
  
  set.seed(123)
  dna_folds <- group_vfold_cv(train_raw, v = 10, group = dim_patient_id)
  
  
  dna_folds %>%
    mutate(
      # 1. Check the Validation (Holdout) Folds
      val_prop = map_dbl(splits, function(split) {
        df <- assessment(split)
        mean(df$dna_outcome == "DNA") # Adjust "missed" to your positive class label
      }),
      
      # 2. Check the Training Folds
      train_prop = map_dbl(splits, function(split) {
        df <- analysis(split)
        mean(df$dna_outcome == "DNA")
      })
    ) %>%
    # 3. Keep only the fold ID and the calculated proportions
    select(id, train_prop, val_prop)
  
  #cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK", master = "localhost")
  #doParallel::registerDoParallel(cl)
  
  registerDoFuture()
  plan(multisession, workers = num_workers, maxSizeOfObjects = 2000 * 1024^2)
  
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
  

  
  # best_params
  best_params <- fits %>%
    select_best(metric = "pr_auc")
  
  
  
  hyper_param_trace <- fits %>%
    collect_metrics() %>%
    filter(.metric == "pr_auc") %>%
    select(mean, mtry:min_n) %>%
    pivot_longer(mtry:min_n, values_to = "value", names_to = "parameter") %>%
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(alpha = 0.8, show.legend = FALSE) +
    facet_wrap( ~ parameter, scales = "free_x") +
    labs(x = NULL, y = "PR AUC");hyper_param_trace
  
  para_coord_plot <- fits %>%
    collect_metrics() %>%
    filter(.metric == "pr_auc") %>%
    select(mtry, trees, min_n, mean) %>%
    ggparcoord(
      columns = 1:3,
      # The columns for your hyperparameters
      groupColumn = 4,
      # Color the lines by the 'mean' (PR-AUC) column
      scale = "uniminmax",
      # Scales each column 0-1 so they are visually comparable
      showPoints = TRUE,
      # Adds dots on the axes for each model
      alphaLines = 0.7         # Makes lines slightly transparent to see overlaps
    ) +
    scale_color_viridis_c(option = "viridis", name = "Mean PR-AUC") +
    theme_minimal() +
    labs(
      title = "Parallel Coordinates Plot of Random Forest Tuning",
      subtitle = "Higher PR-AUC paths are highlighted in yellow/green",
      x = "Hyperparameters",
      y = "Normalized Scale (0 to 1)"
    ) +
    theme(
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5),
      legend.position = "right"
    );para_coord_plot
  
  
  # Performance with best params
  
  roc_plot <- collect_predictions(fits, parameters = best_params) |> (\(preds) {
    auc_val <- preds |>
      roc_auc(truth = dna_outcome, .pred_DNA) |>
      pull(.estimate)
    
    preds |>
      roc_curve(truth = dna_outcome, .pred_DNA) |>
      ggplot(aes(x = 1 - specificity, y = sensitivity)) +
      geom_abline(slope = 1,
                  linetype = 2,
                  alpha = 0.4) +
      geom_path(linewidth = 1, color = "midnightblue") +
      coord_equal() +
      theme_bw() +
      labs(
        title = "ROC curve",
        subtitle = paste0("ROC AUC: ", round(auc_val, 3)),
        x = "False positive rate (1 - Specificity)",
        y = "True positive rate (Sensitivity)"
      )
  })();roc_plot
  
  pr_plot <- collect_predictions(fits, parameters = best_params) |> (\(preds) {
    pr_val <- preds |>
      pr_auc(truth = dna_outcome, .pred_DNA) |>
      pull(.estimate)
    
    preds |>
      pr_curve(truth = dna_outcome, .pred_DNA) %>%
      ggplot(aes(x = recall, y = precision)) +
      geom_path(linewidth = 1, color = "midnightblue") +
      geom_hline(yintercept = 0.05,
                 lty = 2,
                 color = "red") + # Baseline
      coord_equal() +
      theme_bw() +
      labs(
        title = "PR Curve",
        subtitle = paste0("PR AUC: ", round(pr_val, 3)),
        x = "Recall (proportion of DNA identified)",
        y = "Precision (reliability of the prediction)"
      )
  })();pr_plot
  
  
  list(
    tune_res = fits,
    best_params = best_params,
    param_trace = hyper_param_trace,
    para_coord_plot = para_coord_plot,
    pr_plot = pr_plot,
    roc_plot = roc_plot
  )
})

saveRDS(model_tune, "data/rf_tune.RDS")
saveRDS(fits, "S:/Finance/Shared Area/BNSSG - BI/8 Modelling and Analytics/working/nh/projects/dna_predictor/data/rf_tune.RDS")

