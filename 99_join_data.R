library(tidyverse)

dataset <- readr::read_csv("data/dna_combined_clean.csv")

# filter clinic which only records attendance retroactively 
dataset <- dataset %>%
  filter(clinic_code != "ENTO/ERS")


dataset_kwok <- readr::read_csv("data/data_kwok.csv")


dataset_kwok <- dataset_kwok %>%
mutate(test_train = case_when(
    test_train == "Testing" ~ "Training",
    test_train == "Training" ~ "Testint",
    TRUE ~ test_train # Keeps everything else the same
  ))


setdiff(names(dataset_kwok), names(dataset))
setdiff(names(dataset), names(dataset_kwok))


data_joined <- bind_rows(
dataset_kwok %>%
  select(-any_of(setdiff(names(dataset_kwok), names(dataset)))) %>%
  mutate(src = "k", 
  across(c(local_spec_code, consultation_media_code, gender_code),  as.character),
  across(c(distance_km, index_multiple_deprivation_decile),  as.numeric),
),

dataset %>%
  select(-`...1`) %>%
mutate(src = "n")

)



data_joined$dim_patient_id[is.na(data_joined$dim_patient_id)] <- replicate(sum(is.na(data_joined$dim_patient_id)), uuid::UUIDgenerate())


saveRDS(data_joined, file = "data/data_joined.RDS")





