source("00_libraries.R")
source("01_hyperparams.R")

model_data <- local({

  ds <- read.csv("data/dna_data_24042026.csv")

  # test train split
  idx <- sample(seq_len(nrow(ds)), size = floor(nrow(ds) * prop_tt))

  ds$dna_outcome <- ifelse(ds$AttendedStatusCode == "3", 1, 0)
  ds$imd <- as.character(ds$Index_Multiple_Deprivation_Decile)
  ds$imd <- ifelse(is.na(ds$imd), "unknown", ds$imd)

  target <- "dna_outcome"

  vars <- c(
    "Local_Spec_code",
    "National_Spec_Code",
    "AppointmentType",
    "distance_km",
    "NFA_Ind",
    "Age.Group",
    "AgeAtAppointment",
    "Ethnicity",
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
    "SITE_CODE",
    "prev_dna_LY",
    "Appt_Hour",
    "Appt_DOW",
    "Appt_Month",
    "Appt_Wknd_Ind",
    "Referral_Urgency",
    "Lead_Time_days",
    "CLINIC_CODE",
    "CLINIC_LOCATION",
    "imd"
  )

  features <- ds %>%
    select(any_of(vars)) %>%
    rename_with(str_to_lower) %>%
    mutate(across(where(is.character), factor))

  # Lump and count neurodiv
  features <- features %>%
    mutate(
      neurodiv = if_any(starts_with("a_"), ~ .x == 1),
      neurodiv_count = rowSums(select(., starts_with("a_")))
    )

  # Temporal features
  features <- features %>%
    mutate(
      appt_month = factor(appt_month %>% substring(1, 10) %>% as.Date(format = "%d/%m/%Y") %>% format("%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      is_school_run = ifelse(appt_hour %in% c(8, 9, 15, 16), 1, 0),
      has_dna_history = prev_dna_ly > 0
    )

  # Spatial features

  features <- features %>%
    mutate(
      is_local = distance_km < 5,
      is_very_far = distance_km > 20
    )
  
  # LLM encoding of high cardinality vars
  llm_train <- bind_cols(target = ds[idx, target], features[idx, ])
  llm_test <- bind_cols(target = ds[-idx, target], features[-idx, ])

  train_lmm_encoders <- function(data, target_col, high_card_cols) {
    map(high_card_cols, function(col) {
      form <- as.formula(paste0(target_col, " ~ (1 | ", col, ")"))
      fit <- glmer(form, data = data, family = binomial)

      # Extract effects and global intercept
      effects <- ranef(fit)[[1]]
      colnames(effects) <- "rel_adj"
      effects[[col]] <- rownames(effects)

      list(
        mapping = effects,
        global_intercept = fixef(fit)[[1]],
        col_name = col
      )
    }) %>%
      set_names(high_card_cols)
  }

  apply_lmm <- function(x, col_name, encoder_list) {
    enc <- encoder_list[[col_name]]

    res <- enc$mapping$rel_adj[match(x, enc$mapping[[col_name]])]

    enc$global_intercept + ifelse(is.na(res), 0, res)
  }

  high_card_cols <- features %>%
    select(where(is.factor)) %>%
    summarise(across(everything(), n_distinct)) %>%
    pivot_longer(cols = everything()) %>%
    filter(value > 20) %>%
    pull(name)

  # LLM encoding on using training data
  my_encoders <- train_lmm_encoders(llm_train, "target", high_card_cols)

  model_matrix <- bind_cols(ds[target], features)
  # apply encoding
  model_matrix <- model_matrix %>%
    mutate(across(all_of(high_card_cols), 
  ~apply_lmm(.x, cur_column(), my_encoders),
  .names = "{.col}_encoded")) %>%
    select(!any_of(high_card_cols))
  
  # imputation
  model_matrix <- model_matrix %>%
    mutate(dna_outcome = factor(dna_outcome)) %>%
    mutate(across(where(is.factor), \(x) coalesce(x, "missing"))) %>%
    mutate(across(where(is.numeric), \(x) coalesce(x, mean(x, na.rm = TRUE))))
  
  ds_train <- model_matrix[idx, ]
  ds_test <- model_matrix[-idx, ]

  list(train = ds_train, test = ds_test)
})
