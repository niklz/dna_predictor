# -------------------------------------------------------------------------
# 1. Setup & Configuration
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()

# -------------------------------------------------------------------------
# 2. Load & Engineer Base Data
# -------------------------------------------------------------------------
dataset <- readRDS("data/data_joined.RDS") 

model_data <- dataset %>%
  mutate(
    # Create target outcome
    dna_outcome = factor(
      ifelse(attended_status_code == "3", "DNA", "attended"),
      levels = c("DNA", "attended") # Ensure 'DNA' is the first level if it's the event of interest
    ),
    # Handle missing IMD
    imd = coalesce(as.character(index_multiple_deprivation_decile), "unknown")
  ) %>%
  group_by(dim_patient_id) %>%
  # Flag historical missingness
  mutate(patient_profile = if_else(any(dna_outcome == "DNA"), "has_missed", "always_shows")) %>%
  ungroup()

# -------------------------------------------------------------------------
# 3. Split Train / Test (Using Config Vars)
# -------------------------------------------------------------------------
train_raw <- model_data %>%
  filter(test_train == "Training") %>%
  select(all_of(c(conf$target_col, conf$vars, "dim_patient_id")))

test_raw <- model_data %>%
  filter(test_train != "Training") %>%
  select(all_of(c(conf$target_col, conf$vars, "dim_patient_id")))

# -------------------------------------------------------------------------
# 4. Define the Preprocessing Recipe 
# -------------------------------------------------------------------------
train_engineered <- apply_custom_feature_engineering(train_raw)

dna_recipe <- recipe(as.formula(paste(conf$target_col, "~ .")), data = train_engineered) %>%
  
  # 0. Keep ID for tracking predictions later, but hide it from the model
  update_role(dim_patient_id, new_role = "id") %>%
  
  # 1. Handle novel factor levels safely when predicting on future patients
  step_novel(all_nominal_predictors()) %>%
  
  # 2. Handle missing categorical data (excluding 'imd' since it was handled prior)
  step_unknown(all_nominal_predictors(), -imd) %>%
  
  # 3. Lump rare factor levels together dynamically using our config parameter
  step_other(all_nominal_predictors(), threshold = conf$fct_other_prp) %>%
  
  # 4. Target encoding for high-cardinality nominals (GP Practices, Clinics)
  # Notice the dynamic injection of the target column using !!sym()
  step_lencode_mixed(
    any_of(c("clinic_location", "clinic_code", "site_code", "registered_gp_practice")), 
    outcome = vars(!!sym(conf$target_col))
  ) %>%
  
  # 5. Remove zero and near-zero variance predictors
  step_nzv(all_predictors()) %>%
  
  # 6. Impute any missing numerical data with the median
  step_impute_median(all_numeric_predictors())

# -------------------------------------------------------------------------
# 5. Save Outputs (Base R)
# -------------------------------------------------------------------------
dir.create("data/processed", showWarnings = FALSE)

saveRDS(train_engineered, "data/processed/train_engineered.rds")
saveRDS(test_raw, "data/processed/test_raw.rds")
saveRDS(dna_recipe, "data/processed/dna_recipe.rds")

message("Data prep complete. Train, Test, and Recipe saved to /data/processed/")
