library(tidyverse)
library(tidymodels)
library(embed) # For step_lmer (High-cardinality encoding)
library(themis)
library(lme4)
library(ranger)
library(dplyr)
library(probably)
library(patchwork)

library(future)
library(doFuture)
library(tictoc)

dataset <- readRDS("data/data_joined.RDS")

fct_other_prp <- 0.02

  # --- 1. Data Prep ---
  dataset <- dataset %>%
    mutate(
      dna_outcome = factor(
        ifelse(attended_status_code == "3", "DNA", "attended"),
        levels = c("DNA", "attended")
      ),
      imd = coalesce(as.character(index_multiple_deprivation_decile), "unknown")
    )

  target_col <- "dna_outcome"

  vars <- c(
    "local_spec_code",
    "national_spec_code",
    "appointment_type",
    "distance_km",
    "nfa_ind",
    "age_group",
    "age_at_appointment",
    "ethnicity",
    "a_ld",
    "a_autism",
    "a_interpreter_req_bsl",
    "a_interpreter_req_lang",
    "a_balance",
    "a_cognitive_impairment",
    "a_mobility_restriction",
    "a_hear_vis_impaired",
    "a_dementia",
    "a_depression",
    "a_downs_syndrome",
    "a_long_standing_condition",
    "a_makaton",
    "a_mild_cognitive_impairment",
    "a_memory_impairment",
    "a_mood_disorder",
    "a_other_disability",
    "a_psychosis",
    "a_severe_anxiety",
    "a_wheelchair_user",
    "gender",
    "registered_gp_practice",
    "site_code",
    "prev_dna_ly",
    "appt_hour",
    "appt_dow",
    "appt_month",
    "appt_wknd_ind",
    "referral_urgency",
    "lead_time_days",
    "clinic_code",
    "clinic_location",
    "imd"
  )

  # Separate future/unseen data based on your specific index
  train_raw <- dataset %>%
    filter(test_train == "Training") %>%
    select(all_of(c(target_col, vars)))
  test_raw <- dataset %>%
    filter(test_train != "Training") %>%
    select(all_of(c(target_col, vars)))

  # --- 2. The Recipe (Pre-processing Pipeline) ---
  # Tidymodels handles "knowledge separation" automatically.
  dna_recipe <- recipe(dna_outcome ~ ., data = train_raw) %>%
    step_mutate(
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
    step_rm(appt_hour, lead_time_days, appt_date, appt_month, prev_dna_ly) %>%
    step_novel(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors(), -imd) %>%
    step_other(all_nominal_predictors(), threshold = fct_other_prp) %>%
    step_zv(all_predictors()) %>% 
    step_nzv(all_predictors()) %>% 
    step_lencode_mixed(
      clinic_location,
      clinic_code,
      site_code,
      registered_gp_practice,
      outcome = vars(dna_outcome)
    ) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_downsample(dna_outcome, under_ratio = 1)
  
  dna_recipe_ns <- recipe(dna_outcome ~ ., data = train_raw) %>%
    step_mutate(
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
    step_rm(appt_hour, lead_time_days, appt_date, appt_month, prev_dna_ly) %>%
    step_novel(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors(), -imd) %>%
    step_other(all_nominal_predictors(), threshold = fct_other_prp) %>%
    step_zv(all_predictors()) %>% 
    step_nzv(all_predictors()) %>% 
    step_lencode_mixed(
      clinic_location,
      clinic_code,
      site_code,
      registered_gp_practice,
      outcome = vars(dna_outcome)
    ) %>%
    step_impute_median(all_numeric_predictors()) 

  # --- 3. Model Specification ---
rf_spec <- rand_forest(
    mtry = tune(),
    trees = tune(),
    min_n = tune() 
  ) %>%
    set_engine("ranger",
   num.threads = 4,
   #num.threads = parallel::detectCores(),
    importance = "permutation") %>%
    set_mode("classification")

prepped_features <- prep(dna_recipe) %>% juice() %>% select(-dna_outcome)

rf_grid <- grid_space_filling(
  mtry() %>% finalize(prepped_features),
  min_n(),
  trees(range = c(250, 1000)),
  size = 25
)

set.seed(123)
dna_folds <- vfold_cv(train_raw, v = 10, strata = dna_outcome)

#cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK", master = "localhost")
#doParallel::registerDoParallel(cl)

registerDoFuture()
plan(multisession, workers = 9)

tic("Tidymodels grid tuning")

fits <- tibble(
  id = c("down_sampling", "no_sampling"),
  recipe = list("down_sampling" = dna_recipe, "no_sampling" =dna_recipe_ns)
) %>%
  mutate(
    workflow = map(recipe, \(x) {
      workflow() %>% add_recipe(x) %>% add_model(rf_spec)
    }),
    cv_tune = map(workflow, \(x) {
      tune_grid(
        x,
        resamples = dna_folds,  
        grid = rf_grid,        
        metrics = metric_set(pr_auc, roc_auc),
        control = control_grid(
          save_pred = TRUE,
          save_workflow = TRUE,
          parallel_over = "everything"
        )
      )
    })
  )

toc()

saveRDS(fits, "data/rf_tuning_fits.RDS")


fits %>%
  mutate(best_params = map(cv_tune, \(x) select_best(x, metric = "pr_auc"))) %>% pull(best_params)
  mutate(
    predictions = map2(cv_tune, best_params, \(x, y) {
      x %>% collect_predictions(parameters = y)
    })
  ) %>%
  mutate(
    cal_model = map(predictions, \(x) {
      cal_estimate_logistic(
        x,
        truth = dna_outcome,
        estimate = c(.pred_DNA, .pred_attended),
        event_level = "first"
      )
    })
  ) %>%
  mutate(
    predictions_cal = map2(predictions, cal_model, \(x, y) cal_apply(x, y))
  ) %>%
  mutate(
    pred_prob_dens = map(predictions, \(x) {
      x %>%
        ggplot(aes(x = .pred_DNA)) + # Map fill color to the risk value
        geom_histogram(
          bins = 30, # More granular bins are usually better for probabilities
          color = "white", # Thin white border around bins makes them distinct
          fill = "#9a4b53",
          size = 0.2,
          position = "identity",
          closed = "right" # Controls how boundary cases are handled
        ) +
        scale_x_continuous(
          labels = percent_format(accuracy = 1), # Format x-axis as percentages
          expand = c(0.01, 0), # Reduce empty space on the sides
          breaks = seq(0, 1, by = 0.1) # Force breaks at every 10%
        ) +
        scale_y_continuous(
          labels = comma, # Add commas to high y-axis counts (e.g., 1,000)
          expand = expansion(mult = c(0, 0.1)) # Add 10% space at the top so bars don't touch the edge
        ) +
        # color_palette +
        # --- Theme & Labs ---
        theme_minimal() + # Use a clean, minimal base theme
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "grey40"),
          panel.grid.minor = element_blank(), # Remove minor grid lines for a cleaner look
          legend.position = "none" # Remove the legend (x-axis already explains the colors)
        ) +
        labs(
          title = "Distribution of predicted DNA probability",
          subtitle = "Analysis of patient non-attendance risk scores",
          x = "Predicted DNA probability ",
          y = "Frequency (number of appointments)"
        )
    })
  ) %>%
  mutate(
    pred_prob_dens_cal = map(predictions_cal, \(x) {
      x %>%
        ggplot(aes(x = .pred_DNA)) + # Map fill color to the risk value
        geom_histogram(
          bins = 30, # More granular bins are usually better for probabilities
          color = "white", # Thin white border around bins makes them distinct
          fill = "#9a4b53",
          size = 0.2,
          position = "identity",
          closed = "right" # Controls how boundary cases are handled
        ) +
        scale_x_continuous(
          labels = percent_format(accuracy = 1), # Format x-axis as percentages
          expand = c(0.01, 0), # Reduce empty space on the sides
          breaks = seq(0, 1, by = 0.1) # Force breaks at every 10%
        ) +
        scale_y_continuous(
          labels = comma, # Add commas to high y-axis counts (e.g., 1,000)
          expand = expansion(mult = c(0, 0.1)) # Add 10% space at the top so bars don't touch the edge
        ) +
        # color_palette +
        # --- Theme & Labs ---
        theme_minimal() + # Use a clean, minimal base theme
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "grey40"),
          panel.grid.minor = element_blank(), # Remove minor grid lines for a cleaner look
          legend.position = "none" # Remove the legend (x-axis already explains the colors)
        ) +
        labs(
          title = "Distribution of predicted DNA probability",
          subtitle = "Analysis of patient non-attendance risk scores",
          x = "Predicted DNA probability ",
          y = "Frequency (number of appointments)"
        )
    })
  ) %>%
  mutate(
    roc_auc_cv = map_dbl(predictions, \(x) {
      x %>% roc_auc(truth = dna_outcome, .pred_DNA) %>% pull(.estimate)
    })
  ) %>%
  mutate(
    pr_auc_cv = map_dbl(predictions, \(x) {
      x %>% pr_auc(truth = dna_outcome, .pred_DNA) %>% pull(.estimate)
    })
  ) %>%
  mutate(
    pr_curve = map2(predictions, pr_auc_cv, \(x, y) {
      x %>%
        pr_curve(truth = dna_outcome, .pred_DNA) %>%
        ggplot(aes(x = recall, y = precision)) +
        geom_path(linewidth = 1, color = "midnightblue") +
        geom_hline(yintercept = 0.05, lty = 2, color = "red") + # Baseline
        coord_equal() +
        theme_bw() +
        labs(
          title = "PR Curve",
          subtitle = paste0("PR AUC: ", round(y, 3)),
          x = "Recall (proportion of DNA identified)",
          y = "Precision (reliability of the prediction)"
        )
    })
  ) %>%
  mutate(
    roc_curve = map2(predictions, roc_auc_cv, \(x, y) {
      x %>%
        roc_curve(truth = dna_outcome, .pred_DNA) %>%
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_abline(slope = 1, linetype = 2, alpha = 0.4) +
        geom_path(linewidth = 1, color = "midnightblue") +
        # geom_hline(yintercept = 0.05, lty = 2, color = "red") + # Baseline
        coord_equal() +
        theme_bw() +
        labs(
          title = "ROC curve",
          subtitle = paste0("ROC AUC: ", round(y, 3)),
          x = "False positive rate (1 - Specificity)",
          y = "True positive rate (Sensitivity)"
        )
    })
  ) %>%
  select(id, pred_prob_dens, pred_prob_dens_cal, pr_curve, roc_curve) %>%
  mutate(ptc = pmap(list(pr_curve, roc_curve, pred_prob_dens, pred_prob_dens_cal), \(a, b, c, d) a | b | c | d)) %>%
    pull(ptc) %>%
    wrap_plots(ncol = 1)



  dna_workflow <- workflow() %>%
    add_recipe(dna_recipe) %>%
    add_model(rf_spec)

  # --- 5. Fit & Evaluate ---
  # We need to save predictions to generate the PR Curve later
  cv_results <- fit_resamples(
    dna_workflow,
    resamples = dna_folds,
    metrics = metric_set(pr_auc, roc_auc, precision, recall),
    control = control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

  best_rf <- select_best(cv_results, metric = "pr_auc")

  roc_auc_cv <- cv_results %>%
    collect_predictions(parameters = best_rf) %>%
    roc_auc(truth = dna_outcome, .pred_DNA) %>%
    pull(.estimate)

  pr_auc_cv <- cv_results %>%
    collect_predictions(parameters = best_rf) %>%
    pr_auc(truth = dna_outcome, .pred_DNA) %>%
    pull(.estimate)

cv_results %>%
  collect_predictions(parameters = best_rf) %>%
  pr_curve(truth = dna_outcome, .pred_DNA) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path(linewidth = 1, color = "midnightblue") +
  geom_hline(yintercept = 0.05, lty = 2, color = "red") + # Baseline
  coord_equal() +
  theme_bw() +
  labs(
    title = "PR Curve",
    subtitle = paste0("PR AUC: ", round(pr_auc_cv, 3)),
    x = "Recall (Proportion of No-Shows caught)",
    y = "Precision (Reliability of the prediction)"
  )
  
 cv_results %>%
    collect_predictions(parameters = best_rf) %>%
    roc_curve(truth = dna_outcome, .pred_DNA) %>%
    ggplot(aes(x = 1-specificity, y = sensitivity)) +
    geom_path(linewidth = 1, color = "midnightblue") +
    # geom_hline(yintercept = 0.05, lty = 2, color = "red") + # Baseline
    coord_equal() +
    theme_bw() +
    labs(
      title = "Final Tuned XGBoost PR Curve",
      subtitle = paste0("PR AUC: ", round(pr_auc_cv, 3)),
      x = "False positive rate (1 - Specificity)",
      y = "True positive rate (Sensitivity)"
    )


  # Use training fold predictions to build calibration model
  cv_preds <- cv_results %>%
    collect_predictions(parameters = best_rf)

  cal_model <- cal_estimate_logistic(
    cv_preds,
    truth = dna_outcome,
    estimate = c(.pred_DNA, .pred_attended),
    event_level = "first"
  )

  final_fit <- fit(dna_workflow, data = train_raw)

  dataset %>%
    bind_cols(
      predict(final_fit, new_data = ., type = "prob") %>% 
        cal_apply(cal_model)
    ) %>%
    mutate(
      prediction_rank = percent_rank(.pred_DNA), 
      .by = test_train
    ) %>%
    mutate(
      roc_auc_cv = roc_auc_cv, 
      pr_auc_cv = pr_auc_cv
    )




