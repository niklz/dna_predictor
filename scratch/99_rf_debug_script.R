library(tidymodels)
library(embed) # For step_lmer (High-cardinality encoding)
library(themis)
library(lme4)
library(ranger)
library(dplyr)
library(probably)
library(vip)


dataset <- read.csv("data/dna_combined_clean.csv")

# filter clinic which only records attendance retroactively 
dataset <- dataset %>%
  filter(clinic_code != "ENTO/ERS")



# Tuned values 21/05/2026
min_n <- 27
mtry <- 2
trees <- 843
fct_other_prp <- 0.02

# --- 1. Data Prep ---
dataset <- dataset %>%
  mutate(
    dna_outcome = factor(ifelse(attended_status_code == "3", "DNA", "attended"), 
                         levels = c("DNA", "attended")),
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
train_raw <- dataset %>% filter(test_train == "Training") %>%
  select(all_of(c(target_col, vars)))
test_raw  <- dataset %>% filter(test_train != "Training") %>%
  select(all_of(c(target_col, vars)))

# --- 2. The Recipe (Pre-processing Pipeline) ---
# Tidymodels handles "knowledge separation" automatically.
dna_recipe <- recipe(dna_outcome ~ ., data = train_raw) %>%
    step_mutate(
      appt_date = as.Date(substring(appt_month, 1, 10), format = "%d/%m/%Y"),
      appt_dow = factor(weekdays(appt_date)),
      appt_month_num = as.factor(format(appt_date, "%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      # Use pmax to floor lead_time at 0
      lead_time_days_log = log1p(pmax(0, lead_time_days)),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      appt_hour_sin = sin(2 * pi * appt_hour / 24),
      appt_hour_cos = cos(2 * pi * appt_hour / 24),
      has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
    ) %>%
    step_rm(appt_hour, lead_time_days, appt_date, appt_month, prev_dna_ly) %>%
    # Novel levels catch-all
    step_novel(all_nominal_predictors()) %>%
  step_zv(all_predictors()) %>% # Removes zero-variance predictors
  step_nzv(all_predictors()) %>% # Removes near-zero variance predictors
    # Encoding
    step_lencode_mixed(
      clinic_location,
      clinic_code,
      site_code,
      registered_gp_practice,
      outcome = vars(dna_outcome)
    ) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_unknown(all_nominal_predictors(), -imd) %>%
    step_nzv(all_predictors()) %>%
    step_other(all_nominal_predictors(), threshold = fct_other_prp) %>%
    step_downsample(dna_outcome, under_ratio = 1)



# --- 3. Model Specification ---
rf_spec <- rand_forest(
  mtry = mtry,
  trees = trees,
  min_n = min_n # your 'nodes' hyperparam
) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

set.seed(123)
dna_folds <- vfold_cv(train_raw, v = 10, strata = dna_outcome)

# --- 4. The Workflow (The Container) ---
# This binds the recipe and model so they act as a single unit
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


roc_auc_cv <-  cv_results %>%
  collect_predictions(parameters = best_rf) %>%
  roc_auc(truth = dna_outcome, .pred_DNA) %>%
  pull(.estimate)

pr_auc_cv <-  cv_results %>%
  collect_predictions(parameters = best_rf) %>%
  pr_auc(truth = dna_outcome, .pred_DNA) %>%
  pull(.estimate)
  

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

predictions <- predict(final_fit, dataset, type = "prob")

predictions_calibrated <- cal_apply(predictions, cal_model)
predictions_rank <- percent_rank(predictions_calibrated$.pred_DNA)

bind_cols(mutate(dataset, roc_auc_cv = roc_auc_cv, pr_auc_cv = pr_auc_cv),  predictions_calibrated, prediction_rank = predictions_rank)


vip_plot <- vip(final_fit, num_features = 20, geom = "point") +
  theme_minimal() +
  labs(title = "Top 20 Features Driving Predictions",
       subtitle = "Check top features for potential future-leakage");vip_plot


cv_p <- cv_results %>% 
  collect_predictions(parameters = best_rf) %>% 
  cal_apply(cal_model) %>%
  mutate(source = "Cross-Validation Train") %>% 
  select(.pred_DNA, source)

debug_preds_cal <- predict(final_fit, test_raw, type = "prob") %>%
  cal_apply(cal_model) %>% 
  mutate(source = "Production") %>% 
  select(.pred_DNA, source)

combined_preds <- bind_rows(cv_p, debug_preds_cal)

drift_plot <- ggplot(combined_preds, aes(x = .pred_DNA, fill = source)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Distribution Shift Check: Train vs. Future",
       x = "Predicted Probability of DNA",
       y = "Density");drift_plot


cat("--- Diagnostic 3: Novel Level & Missingness Exposure ---\n")

check_novel_levels <- function(train_col, future_col, name) {
  novel_count <- sum(!(future_col %in% train_col))
  pct_novel <- (novel_count / length(future_col)) * 100
  cat(paste0(name, ": ", novel_count, " unseen levels found in future data (", round(pct_novel, 2), "% of rows)\n"))
}

check_novel_levels(train_raw$clinic_code, test_raw$clinic_code, "clinic_code")
check_novel_levels(train_raw$registered_gp_practice, test_raw$registered_gp_practice, "registered_gp_practice")
check_novel_levels(train_raw$clinic_location, test_raw$clinic_location, "clinic_location")
