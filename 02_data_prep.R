source("00_libraries.R")
source("01_hyperparams.R")

model_data <- local({
dataset <- readRDS("data/data_joined.RDS")
dataset <- dataset %>%
  mutate(
    dna_outcome = factor(
      ifelse(attended_status_code == "3", "DNA", "attended"),
      levels = c("DNA", "attended")
    ),
    imd = coalesce(as.character(index_multiple_deprivation_decile), "unknown")
  ) %>%
  group_by(dim_patient_id) %>%
  # If they missed even one appointment, label them "has_missed", otherwise "always_shows"
  mutate(patient_profile = if_else(any(dna_outcome == "DNA"), "has_missed", "always_shows")) %>%
  ungroup()

dataset
})

# Separate future/unseen data based on your specific index
train_raw <- model_data %>%
  filter(test_train == "Training") %>%
  select(all_of(c(target_col, vars, "dim_patient_id")))
test_raw <- model_data %>%
  filter(test_train != "Training") %>%
  select(all_of(c(target_col, vars, "dim_patient_id")))
