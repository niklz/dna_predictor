library(tidyverse)
library(tidymodels)
library(embed) # For step_lmer (High-cardinality encoding)
library(themis)
library(lme4)
library(ranger)
library(dplyr)
library(probably)
library(patchwork)
library(GGally)

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
    select(all_of(c(target_col, vars, "dim_patient_id")))
  test_raw <- dataset %>%
    filter(test_train != "Training") %>%
    select(all_of(c(target_col, vars, "dim_patient_id")))

  # --- 2. The Recipe (Pre-processing Pipeline) ---
  # Tidymodels handles "knowledge separation" automatically.
  dna_recipe <- recipe(dna_outcome ~ ., data = train_raw) %>%
    update_role(dim_patient_id, new_role = "id") %>%
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
    step_lencode_mixed(any_of(c("clinic_location", "clinic_code", "site_code", "registered_gp_practice")),
      outcome = vars(dna_outcome)
    ) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_downsample(dna_outcome, under_ratio = 1)
  
  dna_recipe_ns <- recipe(dna_outcome ~ ., data = train_raw) %>%
    update_role(dim_patient_id, new_role = "id") %>%
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
    step_lencode_mixed(any_of(c("clinic_location", "clinic_code", "site_code", "registered_gp_practice")),
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
saveRDS(fits, "S:/Finance/Shared Area/BNSSG - BI/8 Modelling and Analytics/working/nh/projects/dna_predictor/data/rf_tuning_fits.RDS")

fits <- readRDS("data/rf_tuning_fits.RDS")

# best_params
fits %>%
  mutate(best_params = map(cv_tune, \(x) select_best(x, metric = "pr_auc"))) %>% pull(best_params)



fit_diag <- fits %>%
  mutate(best_params = map(cv_tune, \(x) select_best(x, metric = "pr_auc"))) %>%
  mutate(param_trace = map(cv_tune, \(x)
    x %>%
  collect_metrics() %>%
  filter(.metric == "pr_auc") %>% 
  select(mean, mtry:min_n) %>%
  pivot_longer(mtry:min_n,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")
  )) %>% 
  mutate(para_coord_plot = map(cv_tune, \(x)
    x %>%
  collect_metrics() %>%
  filter(.metric == "pr_auc") %>%
  select(mtry, trees, min_n, mean) %>%
  ggparcoord(
    columns = 1:3,           # The columns for your hyperparameters
    groupColumn = 4,         # Color the lines by the 'mean' (PR-AUC) column
    scale = "uniminmax",     # Scales each column 0-1 so they are visually comparable
    showPoints = TRUE,       # Adds dots on the axes for each model
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
  )
  )) %>% 
  mutate(
    predictions = map2(cv_tune, best_params, \(x, y) {
      x %>% collect_predictions(parameters = y)
    })
  ) %>% 
    mutate(cal_model = map(predictions, \(x) {
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
  mutate(cal_plot_pre = map(predictions, \(x) cal_plot_breaks(x, truth = dna_outcome,  estimate = c(.pred_DNA, .pred_attended)))) %>% 
  mutate(cal_plot_post = map(predictions_cal, \(x) cal_plot_breaks(x, truth = dna_outcome,  estimate = c(.pred_DNA, .pred_attended)))) %>% 
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
  ) 

fit_diag %>%
  # select(id, pred_prob_dens, pred_prob_dens_cal, pr_curve, roc_curve) %>%
  mutate(ptc = pmap(list(pr_curve, roc_curve, pred_prob_dens, pred_prob_dens_cal), \(a, b, c, d) a | b | c | d)) %>%
    pull(ptc) %>%
    wrap_plots(ncol = 1)


fit_diag %>% mutate(ptc = pmap(list(pr_curve, roc_curve), \(a, b) a | b)) %>%
    pull(ptc) %>%
    wrap_plots(ncol = 1) +
  plot_layout(axes = "collect")
  
fit_diag %>% pull(pred_prob_dens) %>%
      wrap_plots(ncol = 1) +
  plot_layout(axes = "collect")

p1 <- fit_diag %>% pull(cal_plot_pre) %>%
  wrap_plots(ncol = 1) 

p2 <- fit_diag %>% pull(cal_plot_post) %>% 
  wrap_plots(ncol = 1) 

(p1 | p2) 


cv_results <- fit_diag$cv_tune$no_sampling

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
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_abline(slope = 1, linetype = 2, alpha = 0.4) +
        geom_path(linewidth = 1, color = "midnightblue") +
        # geom_hline(yintercept = 0.05, lty = 2, color = "red") + # Baseline
        coord_equal() +
        theme_bw() +
        labs(
          title = "ROC curve",
          subtitle = paste0("ROC AUC: ", round(roc_auc_cv, 3)),
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

  output <- dataset %>%
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



strata_data <- output %>%
  mutate(
    risk_prob = .pred_DNA,
    # Create the buckets
    risk_strata = case_when(
      risk_prob >= quantile(risk_prob, 0.95) ~ "1: High (Top 5%)",
      risk_prob >= quantile(risk_prob, 0.80) ~ "2: Medium (Next 15%)",
      TRUE ~ "3: Low (Bottom 80%)"
    ),
    # Create a 'predicted' class based on the Top 20% being 'DNA'
    # This replaces the missing .pred_class
    custom_pred = factor(ifelse(risk_strata %in% c("1: High (Top 5%)", "2: Medium (Next 15%)"), 
                                "DNA", "attended"),
                         levels = c("DNA", "attended"))
  )



strata_summary <- strata_data %>%
  group_nest(risk_strata) %>%
  mutate(metrics = map(data, ~count(.x, dna_outcome) %>% mutate(prop = n/sum(n)))) %>%
  unnest(metrics)

ggplot(strata_summary, aes(x = risk_strata, y = prop, fill = dna_outcome)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + # Baseline
  theme_minimal() +
  labs(title = "DNA Rate by Risk Strata",
       subtitle = "Red line represents the population average (5%)",
       y = "Proportion of Patients", x = "Risk Tier")


prop <- mean(train_raw$dna_outcome == "DNA")

strata_data %>%
  group_by(risk_strata) %>%
  summarise(
    total_patients = n(),
    actual_dnas = sum(dna_outcome == "DNA"),
    dna_rate = actual_dnas / total_patients
  ) %>%
  ggplot(aes(x = risk_strata, y = dna_rate, fill = risk_strata)) +
  geom_col(show.legend = FALSE) +
  annotate("text", x = 0.5, y = prop, label = "Population Average", 
           vjust = -1, hjust = 0, color = "grey40", fontface = "italic") +
  geom_hline(yintercept = prop, linetype = "dashed", color = "grey40") + # Baseline
  scale_y_continuous(labels = percent_format()) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  theme_minimal() +
  labs(title = "DNA Rate by Risk Strata",
subtitle = paste0("Dashed line represents the overall average rate of ", round(prop * 100, 1), "%"),
       x = "Assigned Risk Tier",
       y = "Actual DNA rate (%)")

ggsave("materials/strata_plot.png",
last_plot(),
bg = "white",
height = 5,
width = 7,
scale = 0.8)




