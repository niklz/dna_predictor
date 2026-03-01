source("00_libraries.R")
source("01_hyperparams.R")

model_data <- local({
  dataset <- read.csv("data/dna_data_full_clean.csv")
  
  ds <- dataset
  
  ds$dna_outcome <- ifelse(ds$AttendedStatusCode == "3", 1, 0)
  # ds$prev_dna <- ifelse(is.na(ds$prev_dna), 0, ds$prev_dna)
  ds$imd <- as.character(ds$Index_Multiple_Deprivation_Decile)
  ds$imd <- ifelse(is.na(ds$imd), "unknown", ds$imd)
  
  target <- "dna_outcome"
  
  vars <- c(
    "Local_Spec_Code",
    "National_Spec_Code",
    "AppointmentType",
    # "ConsultationMedia",
    "distance_km",
    "NFA_Ind",
    "Age.Group",
    "AgeAtAppointment",
    "Ethnicity",
    "imd",
    "a_LD",
    "a_Autism",
    "a_InterpreterReq_BSL",
    "a_InterpreterReq_lang",
    "a_Balance",
    "a_CognitiveImpairment",
    "a_MobilityRestriction",
    "a_HearVisImpaired",
    "a_Dementia",
    "a_Depression",
    "a_DownsSyndrome",
    "a_LongStandingCondition",
    "a_Makaton",
    "a_MildCognitiveImpairment",
    "a_MemoryImpairment",
    "a_MoodDisorder",
    "a_OtherDisability",
    "a_Psychosis",
    "a_SevereAnxiety",
    "a_WheelchairUser",
    "Gender",
    "RegisteredGPPractice",
    "prev_dna_LY",
    "Appt_Hour",
    "Appt_DOW",
    "Appt_Month",
    "Appt_Wknd_Ind",
    "Referral_Urgency",
    "Lead_Time_days",
    "CLINIC_CODE",
    "CLINIC_LOCATION"
  )
  
  
  features <- ds %>%
    select(any_of(vars)) %>%
    rename_with(str_to_lower)

  # Lump neurodiv

  features <- features %>%
    mutate(
    neurodiv = if_any(starts_with("a_"), ~ .x == 1)
  )
  
  # down_sample_factors
  features <- features %>%
    mutate(across(where(is.character), factor)) %>%
    mutate(across(
      where(is.factor),
      \(x) fct_lump_n(x, n = n_lump, other_level = "other")
    ))
  
  # imputation
  features <- features %>%
    mutate(across(where(is.factor), \(x) coalesce(x, "missing"))) %>%
    mutate(across(where(is.numeric), \(x) coalesce(x, mean(x, na.rm = TRUE))))
  
  # # near zero variance filtration
  # features <-features %>%
  #   select(where(~ {
  #     # Count non-NA values
  #     n_total <- sum(!is.na(.x))
      
  #     # If the column is empty or all NAs, drop it
  #     if (n_total == 0) return(FALSE)
      
  #     # Calculate frequency ratio
  #     counts <- table(.x, useNA = "no")
  #     max_freq <- max(counts)
      
  #     (max_freq / n_total) < nzv_freq
  #   }))
  
  model_matrix <- bind_cols(ds[target], features)
  
  # test train split
  
  idx <- sample(seq_len(nrow(model_matrix)), size = floor(nrow(model_matrix) *
                                                            prop_tt))
  ds_train <- model_matrix[idx, ]
  ds_test <- model_matrix[-idx, ]
  
  list(train = ds_train, test = ds_test)
})
