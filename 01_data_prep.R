# -------------------------------------------------------------------------
# 1. Setup & Configuration
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()

local({
  # -------------------------------------------------------------------------
  # 2. Load & Engineer Base Data
  # -------------------------------------------------------------------------
  dataset <- read.csv("data/DNA_20260818.csv")
  
  model_data <- dataset %>%
    mutate(
      # Create target outcome
      dna_outcome = factor(
        ifelse(attended_status_code == "3", "DNA", "attended"),
        levels = c("DNA", "attended") # Ensure 'DNA' is the first level if it's the event of interest
      )
    ) %>%
    group_by(dim_patient_id) %>%
    # Flag historical missingness
    mutate(patient_profile = if_else(any(dna_outcome == "DNA"), "has_missed", "always_shows")) %>%
    ungroup()
  
  # -------------------------------------------------------------------------
  # 3. Split Train / Test (Using Config Vars)
  # -------------------------------------------------------------------------
  # Swap "imd" for the raw database column so it is available for feature engineering
  select_vars <- setdiff(conf$vars, "imd")
  
  train_raw <- model_data %>%
    filter(test_train == "Training") %>%
    select(all_of(c(
      conf$target_col, 
      select_vars, 
      "index_multiple_deprivation_decile", 
      "dim_patient_id"
    )))
  
  test_raw <- model_data %>%
    filter(test_train != "Training") %>%
    select(all_of(c(
      conf$target_col, 
      select_vars, 
      "index_multiple_deprivation_decile", 
      "dim_patient_id"
    )))
  
  # -------------------------------------------------------------------------
  # 4. Define the Preprocessing Recipe
  # -------------------------------------------------------------------------
  # The engineering step now builds 'imd' cleanly on the fly
  train_engineered <- apply_custom_feature_engineering(train_raw)
  
  data_template <- head(train_engineered, 0)
  
  dna_recipe <- build_trial_recipe(
    data_template = data_template,
    target_col    = conf$target_col,
    fct_other_prp = conf$fct_other_prp
  )
  
  # -------------------------------------------------------------------------
  # 5. Extract a Balanced, Stratified 10% Subset for Grid Tuning
  # -------------------------------------------------------------------------
  # With ~1M rows, running 10-fold CV over a grid of 25 is computationally redundant.
  # We slice a representative 20% sample, group-stratified by our target outcome 
  # (dna_outcome) to preserve baseline class prevalence.
  set.seed(42) # Set seed for reproducible sampling
  
  train_engineered_tune <- train_engineered %>%
    group_by(!!sym(conf$target_col)) %>%
    slice_sample(prop = 0.20) %>%
    ungroup()
  
  # -------------------------------------------------------------------------
  # 6. Save Outputs (Base R)
  # -------------------------------------------------------------------------
  dir.create("data/processed", showWarnings = FALSE)
  
  # Save the full 1,000,000-row set for the final model fit
  saveRDS(train_engineered, "data/processed/train_engineered.rds")
  
  # Save the 100,000-row subset strictly for the grid-tuning loop
  saveRDS(train_engineered_tune, "data/processed/train_engineered_tune.rds")
  
  saveRDS(test_raw, "data/processed/test_raw.rds")
  saveRDS(dna_recipe, "data/processed/dna_recipe.rds")
  
  message("Data prep complete. Full train, 20% tune subset, test set, and recipe saved.")
})
