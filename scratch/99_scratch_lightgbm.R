set.seed(123)
library(bonsai)

dna_folds <- vfold_cv(train_raw, v = 10, strata = dna_outcome)


# 1. Model Spec with Tuning Placeholders
lgbm_spec <- boost_tree(
  mtry = tune(),
  trees = 1000,
  min_n = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune()
) %>%
  set_engine("lightgbm", is_unbalance = TRUE) %>% # Handles imbalance internally
  set_mode("classification")

# 2. Creating a Space-Filling Design Grid
# This covers the 'search space' efficiently
lgbm_grid <- grid_space_filling(
  finalize(mtry(), train_raw),
  min_n(),
  tree_depth(),
  learn_rate(),
  loss_reduction(),
  size = 30 # Number of combinations to try
)

lgbm_workflow <- workflow() %>%
  add_recipe(dna_recipe) %>% # Ensure your recipe includes step_downsample(dna_outcome)
  add_model(lgbm_spec)


# Run the tuning
tune_results <- tune_grid(
  lgbm_workflow,
  resamples = dna_folds,
  grid = lgbm_grid,
  metrics = metric_set(pr_auc, roc_auc, sens, precision),
  control = control_grid(save_pred = TRUE)
)

# Pick the best parameters based on Precision-Recall AUC
best_lgbm_params <- tune_results %>%
  select_best(metric = "pr_auc")

# Finalize the workflow
final_lgbm_wf <- lgbm_workflow %>%
  finalize_workflow(best_lgbm_params)


# 2. Run the workflow across all folds
# Tidymodels is smart: it will downsample inside each fold's training set, 
# but evaluate on that fold's UNTOUCHED (imbalanced) validation set.
cv_results <- fit_resamples(
  final_lgbm_wf,
  resamples = dna_folds,
  metrics = metric_set(pr_auc, roc_auc, precision, recall),
  control = control_resamples(save_pred = TRUE) # Crucial for PR curves
)



# 3. Collect the "Real World" PR Curve data
cv_results %>%
  collect_predictions() %>%
  pr_curve(truth = dna_outcome, .pred_DNA) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  theme_minimal() +
  labs(title = "Cross-Validated PR Curve (Real-World Prevalence)")


cv_results %>%
  collect_predictions() %>%
  pr_auc(truth = dna_outcome, .pred_DNA) 

