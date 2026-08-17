library(tidymodels)
library(embed) # For step_lmer (High-cardinality encoding)
library(themis)
library(lme4)
library(ranger)
library(dplyr)

dataset <- read.csv("data/dna_combined_clean.csv")

dataset <- local({
  nodes <- 100
  mtry <- 7
  fct_other_prp <- 0.02

  # --- 1. Data Prep ---
  dataset <- dataset %>%
    mutate(
      dna_outcome = factor(
        ifelse(attended_status_code == "3", "DNA", "attended"),
        levels = c("DNA", "attended")
      ),
      imd = coalesce(as.character(index_multiple_deprivation_decile), "unknown"),
      imd = factor(imd, levels = as.character(1:10))
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
      # Use pmax to floor lead_time at 0
      lead_time_days_log = log1p(pmax(0, lead_time_days)),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      appt_hour_sin = sin(2 * pi * appt_hour / 24),
      appt_hour_cos = cos(2 * pi * appt_hour / 24),
      has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
    ) %>%
    # step_discretize(
    #   lead_time_days,
    #   num_breaks = 4,
    #   min_unique = 5) %>%
    step_rm(appt_hour, lead_time_days, appt_date, appt_month) %>%
    # Novel levels catch-all
    step_novel(all_nominal_predictors()) %>%
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
    trees = 500,
    min_n = nodes # your 'nodes' hyperparam
  ) %>%
    set_engine("ranger", importance = "permutation") %>%
    set_mode("classification")

  # --- 4. The Workflow (The Container) ---
  # This binds the recipe and model so they act as a single unit
  dna_workflow <- workflow() %>%
    add_recipe(dna_recipe) %>%
    add_model(rf_spec)

  # --- 5. Fit & Evaluate ---
  # This fits the recipe (calculates means/encoders) and the model ONLY on train_raw
  final_fit <- fit(dna_workflow, data = train_raw)

# To get the AUC, we need the outcome variable AFTER downsampling
# 'bake' the recipe on the training data to see what actually went into the model
train_processed <- bake(prep(dna_recipe), new_data = NULL)

# Get OOB performance
oob_probs <- final_fit %>%
  extract_fit_engine() %>%
  pluck("predictions")

# Now 'truth' and 'estimate' will have the same length!
auc_oob <- roc_auc_vec(
  truth = train_processed$dna_outcome,
  estimate = oob_probs[, "DNA"]
)

# Return the full dataset with predictions attached
bind_cols(
  mutate(dataset, auc_oob = auc_oob),
  predict(final_fit, dataset, type = "prob")
)
})

output <- dataset




