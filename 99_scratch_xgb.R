dataset <- read.csv("data/dna_combined_clean.csv")

# filter clinic which only records attendance retroactively 
dataset <- dataset %>%
  filter(clinic_code != "ENTO/ERS")


library(xgboost)
library(probably)
library(tidymodels)
library(embed)
library(themis)

fct_other_prp <- 0.02

# --- 1. Data Prep ---
dataset <- dataset %>%
  mutate(
    dna_outcome = factor(ifelse(attended_status_code == "3", "DNA", "attended"), 
                         levels = c("DNA", "attended")),
    imd = coalesce(as.character(index_multiple_deprivation_decile), "unknown"),
    imd = factor(imd, levels = c(as.character(1:10), "unknown"))
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
xgb_rec <- recipe(dna_outcome ~ ., data = train_raw) %>%
  step_mutate(
    appt_date = as.Date(substring(appt_month, 1, 10), format = "%d/%m/%Y"),
    appt_month_num = as.factor(format(appt_date, "%m")), 
    lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
    # Use pmax to floor lead_time at 0
    lead_time_days_log = log1p(pmax(0, lead_time_days)),
    is_morning = ifelse(appt_hour < 12, 1, 0),
    appt_hour_sin = sin(2 * pi * appt_hour / 24),
    appt_hour_cos = cos(2 * pi * appt_hour / 24),
    has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
  )%>%
  step_rm(appt_hour, lead_time_days, appt_date, appt_month, prev_dna_ly) %>%
  step_novel(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors(), -imd) %>%
  # 1. Encode high-cardinality first (like you already were)
  step_lencode_mixed(
    clinic_location, clinic_code, site_code, registered_gp_practice,
    outcome = vars(dna_outcome)
  ) %>%
  # 2. Convert ALL remaining factors to 0/1 dummy variables (CRITICAL FOR XGBOOST)
  step_dummy(all_nominal_predictors()) %>% 
  # 3. Clean up and downsample
  step_nzv(all_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>% 
  step_downsample(dna_outcome, under_ratio = 1)


# --- 3. Model Specification ---
xgb_spec <- boost_tree(
  trees = 1000,               # Keep this fixed to save time
  tree_depth = tune(),       # How deep the trees go
  min_n = tune(),            # Minimum data in a leaf
  loss_reduction = tune(),   # "Gamma" - prevents overfitting
  sample_size = tune(),      # Stochasticity (subsample)
  mtry = tune(),             # Number of columns per tree
  learn_rate = tune()        # Step size shrinkage
) %>%
  set_engine("xgboost", nthread = parallel::detectCores() - 1, tree_method = "hist") %>%
  set_mode("classification")

xgb_workflow <- workflow() %>%
  add_recipe(xgb_rec) %>% # Uses your existing recipe with downsampling
  add_model(xgb_spec)

# Create a Latin Hypercube grid
set.seed(789)
xgb_grid <- grid_space_filling(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(range = c(5, 15)), train_raw),
  learn_rate(range = c(-2, -0.5)),
  size = 25 
)

dna_folds <- vfold_cv(train_raw, v = 5, strata = dna_outcome)

# Run tuning
xgb_res <- tune_grid(
  xgb_workflow,
  resamples = dna_folds,
  grid = xgb_grid,
  metrics = metric_set(pr_auc, roc_auc),
  control = control_grid(save_pred = TRUE, parallel_over = "everything")
)

saveRDS(xgb_res, "data/xgb_res.RDS")

show_best(xgb_res, metric = "pr_auc")

# Select the winner and finalize
best_xgb <- select_best(xgb_res, metric = "pr_auc")
final_xgb_wf <- finalize_workflow(xgb_workflow, best_xgb)

# Fit one last time on all training data
xgb_fit <- fit(final_xgb_wf, data = train_raw)



roc_auc_cv <-  xgb_res %>%
  collect_predictions() %>%
  roc_auc(truth = dna_outcome, .pred_DNA) %>%
  pull(.estimate)

pr_auc_cv <-  xgb_res %>%
  collect_predictions() %>%
  pr_auc(truth = dna_outcome, .pred_DNA) %>%
  pull(.estimate)


bind_cols(mutate(dataset, roc_auc_cv = roc_auc_cv, pr_auc_cv = pr_auc_cv),  predict(xgb_fit, dataset, type = "prob"))




xgb_res %>%
  collect_predictions() %>%
  pr_curve(truth = dna_outcome, .pred_DNA) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path(linewidth = 1, color = "midnightblue") +
  geom_hline(yintercept = 0.05, lty = 2, color = "red") + # Baseline
  coord_equal() +
  theme_bw() +
  labs(
    title = "Final Tuned XGBoost PR Curve",
    subtitle = paste0("PR AUC: ", round(pr_auc_cv, 3)),
    x = "Recall (Proportion of No-Shows caught)",
    y = "Precision (Reliability of the prediction)"
  )



  xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "pr_auc") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")


prop <- mean(train_raw$dna_outcome == "DNA")

xgb_res %>%
  collect_predictions() %>%
  pr_curve(truth = dna_outcome, .pred_DNA) %>%
  mutate(
    precG = (precision - prop) / ((1 - prop) * precision),
    recG  = (recall - prop) / ((1 - prop) * recall)
  ) %>%
  filter(is.finite(precG), is.finite(recG)) %>%
  # We usually only care about the part of the curve where we beat a random guess
  filter(precG >= 0, recG >= 0) %>%
  ggplot(aes(x = recG, y = precG)) +
  geom_path(linewidth = 1, color = "midnightblue") +
  coord_equal() +
  theme_bw() +
  labs(
    title = "Final Tuned XGBoost PR Curve",
    # subtitle = paste0("PR AUC: ", round(pr_auc_cv, 3)),
    x = "Recall Gain (Proportion of No-Shows caught)",
    y = "Precision Gain (Reliability of the prediction)"
  )


xgb_res %>%
  collect_predictions() %>%
  gain_curve(truth = dna_outcome, .pred_DNA) %>%
  autoplot()



# RF comparison

rf_rec <- recipe(dna_outcome ~ ., data = train_raw) %>%
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
    step_rm(appt_hour, lead_time_days, appt_date, appt_month, prev_dna_ly) %>%
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



  rf_spec <- rand_forest(
    mtry = tune(),
    trees = tune(),
    min_n = tune() 
  ) %>%
    set_engine("ranger", importance = "permutation") %>%
    set_mode("classification")

rf_grid <- grid_space_filling(
  mtry() %>% finalize(train_raw %>% select(-dna_outcome)),
  min_n(),
  trees(range = c(250, 1000)),
  size = 25 
)

 rf_workflow <- workflow() %>%
    add_recipe(rf_rec) %>%
    add_model(rf_spec)


doParallel::registerDoParallel()
rf_res <- tune_grid(
  rf_workflow,
  resamples = dna_folds,
  grid = rf_grid,
  metrics = metric_set(pr_auc, roc_auc),
  control = control_grid(save_pred = TRUE, parallel_over = "everything")
)

saveRDS(rf_res, "data/rf_res_2.RDS")


# Select the winner and finalize
best_rf <- select_best(rf_res, metric = "pr_auc")
final_rf_wf <- finalize_workflow(rf_workflow, best_rf)

# Fit one last time on all training data
rf_fit <- fit(final_rf_wf, data = train_raw)

rf_fit %>% vip::vip()

# Outputs

rf_res <- readRDS("data/rf_res.RDS")
best_rf <- select_best(rf_res, metric = "pr_auc")
xgb_res <- readRDS("data/xgb_res.RDS")
best_xgb <- select_best(xgb_res, metric = "pr_auc")

rf_res %>%
  collect_metrics() %>%
  filter(.metric == "pr_auc") %>%
  select(mean, mtry: min_n) %>%
  pivot_longer(mtry: min_n,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "PR AUC")


bind_rows(
  rf_res %>%
  collect_predictions(parameters = best_rf) %>%

    mutate(model = "Random Forest"),
  xgb_res %>%
  collect_predictions(parameters = best_xgb) %>%
    mutate(model = "XGBoost")
) %>%
  group_by(model) %>%
  gain_curve(truth = dna_outcome, .pred_DNA) %>%
  autoplot() +
  paletteer::scale_color_paletteer_d("nationalparkcolors::Acadia", direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(colour = "",
       title = "Gain Curve for DNA model",
  subtitle = "How quickly the models identify DNAs compared to random selection",
  x = "Percentage of patients contacted",
  y = "Percentage of total DNAs captured",
  caption = "A steeper curve indicates a more effective risk stratification model.")


ggsave("materials/gain_plot.png",
last_plot(),
bg = "white",
height = 8,
width = 7,
scale = 0.8)





bind_rows(
  rf_res %>%
  collect_predictions(parameters = best_rf) %>%

    mutate(model = "Random Forest"),
  xgb_res %>%
  collect_predictions(parameters = best_xgb) %>%
    mutate(model = "XGBoost")
) %>%
  group_by(model) %>%
  pr_curve(truth = dna_outcome, .pred_DNA) %>%
  autoplot()

bind_rows(
  rf_res %>%
  collect_predictions(parameters = best_rf) %>%

    mutate(model = "Random Forest"),
  xgb_res %>%
  collect_predictions(parameters = best_xgb) %>%
    mutate(model = "XGBoost")
) %>%
  group_by(model) %>%
  roc_curve(truth = dna_outcome, .pred_DNA) %>%
  autoplot()

rf_res %>%
  collect_predictions(parameters = best_rf) %>%
  # Create 10 equal-sized groups based on model probability
  mutate(bin = ntile(.pred_DNA, 5)) %>% 
  group_by(bin) %>%
  summarize(
    total_appointments = n(),
    actual_dnas = sum(dna_outcome == "DNA"),
    precision = actual_dnas / total_appointments,
    one_in_x = 1 / precision
  ) %>%
  arrange(desc(bin))

rf_res %>%
  collect_predictions(parameters = best_rf) %>%
  mutate(
    # Define your breaks based on percentiles of the predicted probabilities
    percentile_bin = cut(
      .pred_DNA,
      breaks = quantile(.pred_DNA, probs = c(0, 0.5, 0.85, 0.9, 0.95, 0.99, 1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("0-50%", "50-85%", "85-90%", "90-95%", "95-99%", "99-100%")
    )
  )%>%
  group_by(percentile_bin) %>%
  summarize(
    total = n(),
    actual_dnas = sum(dna_outcome == "DNA"),
    actual_attended = total - actual_dnas,
    precision = actual_dnas / total,
    one_in_x = 1 / precision
  ) %>%
  # Reshape for ggplot stacking
  tidyr::pivot_longer(cols = c(actual_dnas, actual_attended), 
               names_to = "outcome_type", 
               values_to = "count") %>%
  {
  ggplot(., aes(x = percentile_bin, y = count, fill = outcome_type)) +
  # Use position = "fill" for proportional/percentage stacking
  geom_col(position = "fill", alpha = 0.8) +
  # Annotation for "1 in X" at the very top
  geom_text(data = . %>% filter(outcome_type == "actual_dnas"),
            aes(y = 1.05, label = paste0("1 in ", round(one_in_x, 0))),
            size = 3.5, fontface = "bold", color = "grey20") +
  # Labels for counts inside the bars
  geom_text(aes(label = comma(count)), 
            position = position_fill(vjust = 0.5), 
            color = "white", size = 3) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(
    values = c("actual_attended" = "#2a679d", "actual_dnas" = "#E76F51"),
    labels = c("Attended", "DNA (Did Not Attend)")
  ) +
  theme_minimal() +
  labs(
    title = "Proportional DNA outcomes by risk percentile",
    subtitle = "Labels indicate 1-in-X frequency; bar counts show raw volume",
    x = "Predicted risk percentile (Model Score)",
    y = "Proportion of appointments",
    fill = "Outcome"
  ) +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom")
  }

ggsave("materials/risk_percentile_plot.png",
last_plot(),
bg = "white",
height = 8,
width = 7,
scale = 0.8)


roc_auc_cv <-  rf_res %>%
  collect_predictions(parameters = best_rf) %>%
  roc_auc(truth = dna_outcome, .pred_DNA) %>%
  pull(.estimate)

pr_auc_cv <-  rf_res %>%
  collect_predictions(parameters = best_rf) %>%
  pr_auc(truth = dna_outcome, .pred_DNA) %>%
  pull(.estimate)


output <- bind_cols(mutate(dataset, roc_auc_cv = roc_auc_cv, pr_auc_cv = pr_auc_cv),  predict(rf_fit, dataset, type = "prob"))

threshold_data <- rf_res %>%
  collect_predictions(parameters = best_rf) %>%
  threshold_perf(dna_outcome, .pred_DNA, thresholds = seq(0, 1, by = 0.0025)) |>
  filter(.metric != "distance") |>
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold <- threshold_data |>
  filter(.metric == "j_index") |>
  filter(.estimate == max(.estimate)) |>
  pull(.threshold)

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric, alpha = group)) +
  geom_line() +
  theme_minimal() +
  scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30") +
  labs(
    x = "'Good' Threshold\n(above this value is considered 'good')",
    y = "Metric Estimate",
    title = "Balancing performance by varying the threshold",
    subtitle = "Sensitivity or specificity alone might not be enough!\nVertical line = Max J-Index"
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


# Confusion Matrix per Strata
# strata_data %>%
#   group_by(risk_strata) %>%
#   conf_mat(truth = dna_outcome, estimate = custom_pred) %>%
#   mutate(plot = map(conf_mat, autoplot, type = "heatmap")) %>% 
#   autoplot(type = "heatmap") +
#   facet_wrap(~risk_strata) +
#   labs(title = "Confusion Matrix Heatmap by Risk Strata")



output %>% 
  ggplot(aes(x = .pred_DNA)) + # Map fill color to the risk value
  geom_histogram(
    bins = 30,                # More granular bins are usually better for probabilities
    color = "white",          # Thin white border around bins makes them distinct
    fill = "#9a4b53",
    size = 0.2,               
    position = "identity",    
    closed = "right"          # Controls how boundary cases are handled
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1), # Format x-axis as percentages
    expand = c(0.01, 0),                  # Reduce empty space on the sides
    breaks = seq(0, 1, by = 0.1)          # Force breaks at every 10%
  ) +
  scale_y_continuous(
    labels = comma,                       # Add commas to high y-axis counts (e.g., 1,000)
    expand = expansion(mult = c(0, 0.1)) # Add 10% space at the top so bars don't touch the edge
  ) +
  # color_palette +
  # --- Theme & Labs ---
  theme_minimal() + # Use a clean, minimal base theme
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank(), # Remove minor grid lines for a cleaner look
    legend.position = "none"            # Remove the legend (x-axis already explains the colors)
  ) +
  labs(
    title = "Distribution of predicted DNA probability",
    subtitle = "Analysis of patient non-attendance risk scores",
    x = "Predicted DNA probability ",
    y = "Frequency (number of appointments)"
  )

ggsave("materials/risk_hist_plot.png",
last_plot(),
bg = "white",
height = 5,
width = 7,
scale = 0.8)



output %>% 
  ggplot(aes(x = .pred_DNA, fill = .pred_DNA < 0.18)) + # Map fill color to the risk value
  geom_histogram(
    bins = 30,                # More granular bins are usually better for probabilities
    color = "white",          # Thin white border around bins makes them distinct
    # fill = "#9a4b53",
    size = 0.2,               
    position = "identity",    
    closed = "right"          # Controls how boundary cases are handled
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1), # Format x-axis as percentages
    expand = c(0.01, 0),                  # Reduce empty space on the sides
    breaks = seq(0, 1, by = 0.1)          # Force breaks at every 10%
  ) +
  scale_y_continuous(
    labels = comma,                       # Add commas to high y-axis counts (e.g., 1,000)
    expand = expansion(mult = c(0, 0.1)) # Add 10% space at the top so bars don't touch the edge
  ) +
  # color_palette +
  # --- Theme & Labs ---
  theme_minimal() + # Use a clean, minimal base theme
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank(), # Remove minor grid lines for a cleaner look
    legend.position = "none"            # Remove the legend (x-axis already explains the colors)
  ) +
  labs(
    title = "Distribution of predicted DNA probability",
    subtitle = "Analysis of patient non-attendance risk scores",
    x = "Predicted DNA probability ",
    y = "Frequency (number of appointments)"
  )




# ENTO/ERS is driving the low-liklihood
output %>%
  mutate(flag = ifelse(.pred_DNA < 0.18, "low", "normal")) %>%
  select(-c(appt_dttm, dim_patient_id)) %>%
  gtsummary::tbl_summary(by = flag, include = c(clinic_code))




output_2 <- rf_rec %>%
  prep() %>%
  bake(new_data = dataset) %>%
  mutate(prev_dna_ly = dataset$prev_dna_ly, id = dataset$dim_patient_id) %>%
  bind_cols(predict(rf_fit, dataset, type = "prob"))

glimpse(output_2)


output_2 %>%
  filter(.pred_DNA > 0.18) %>%
  mutate(
    lead_time = expm1(lead_time_days_log),
    liklihood_cut = factor(
      cut(
        .pred_DNA,
        quantile(.pred_DNA, c(0, 0.5, .9, .99, 1)),
        include.lowest = TRUE
      ),
      labels = c(
        "Bottom 50%",
        "50-90%",
        "90-99%",
        "99%+"
      )
    )
  ) %>%
  # mutate(liklihood_qtile = ntile(.pred_DNA, 10)) %>%
  gtsummary::tbl_summary(
    by = liklihood_cut,
    include = c(
      age_group,
      has_dna_history,
      prev_dna_ly,
      imd,
      appt_dow,
      is_morning,
      distance_km,
      lead_over_30,
      lead_time
    ),
    label = list(
          age_group = "Age",
          has_dna_history = "Any DNA history",
          prev_dna_ly = "Number DNA (last year)",
          imd = "IMD",
          appt_dow = "Day of week",
          is_morning = "Morning",
          lead_over_30 = "+30 day lead-time",
          lead_time = "Lead time",
          distance_km = "Distance (km)"
        )
  )


output_2 %>%
  filter(.pred_DNA > 0.18) %>%
  mutate(
    lead_time = expm1(lead_time_days_log),
    liklihood_cut = factor(
      cut(
        .pred_DNA,
        quantile(.pred_DNA, c(0, 0.5, .9, .99, 1)),
        include.lowest = TRUE
      ),
      labels = c(
        "Bottom 50%",
        "50-90%",
        "90-99%",
        "99%+"
      )
    )
  ) %>%
  distinct(id, liklihood_cut) %>%
  group_by(group = liklihood_cut) %>%
  sample_n(50) %>%
  select(-liklihood_cut) %>%
  show_in_excel()
