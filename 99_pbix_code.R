dataset <- read.csv("data/dna_combined.csv")

library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(stringr)
library(forcats)
library(ranger)


## hyper params for modelling

# number of levels to lump factors into
n_lump <- 10
# test/train split
prop_tt <- 0.7
# quick feature importance rows to train
n_featimp <- 5E3
# near zero variance frequency
nzv_freq <- 0.80


# RF hyperparams

mtry <- 10
nodes <- 10

# data pre-processing and feature engineering

dataset$dna_outcome <- ifelse(dataset$AttendedStatusCode == "3", 1, 0)
dataset$imd <- as.character(dataset$Index_Multiple_Deprivation_Decile)
dataset$imd <- ifelse(is.na(dataset$imd), "unknown", dataset$imd)

# training indicies
idx <- which(dataset$test_train == "Training")

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

features <- dataset %>%
  select(any_of(vars)) %>%
  rename_with(str_to_lower) %>%
  mutate(across(where(is.character), factor))

# Lump and count neurodiv
features <- features %>%
  mutate(neurodiv = if_any(starts_with("a_"), ~ .x == 1),
         neurodiv_count = rowSums(select(., starts_with("a_"))))

# Temporal features
features <- features %>%
  mutate(
    appt_month = factor(
      appt_month %>% substring(1, 10) %>% as.Date(format = "%d/%m/%Y") %>% format("%m")
    ),
    lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
    is_morning = ifelse(appt_hour < 12, 1, 0),
    is_school_run = ifelse(appt_hour %in% c(8, 9, 15, 16), 1, 0),
    has_dna_history = prev_dna_ly > 0
  )

# Spatial features

features <- features %>%
  mutate(is_local = distance_km < 5, is_very_far = distance_km > 20)

# LLM encoding of high cardinality vars
llm_train <- bind_cols(target = dataset[idx, target], features[idx, ])
llm_test <- bind_cols(target = dataset[-idx, target], features[-idx, ])

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

model_matrix <- bind_cols(dataset[target], features)
# apply encoding
model_matrix <- model_matrix %>%
  mutate(across(
    all_of(high_card_cols),
    ~ apply_lmm(.x, cur_column(), my_encoders),
    .names = "{.col}_encoded"
  )) %>%
  select(!any_of(high_card_cols))

# imputation
model_matrix <- model_matrix %>%
  mutate(dna_outcome = factor(dna_outcome)) %>%
  mutate(across(where(is.factor), \(x) coalesce(x, "missing"))) %>%
  mutate(across(where(is.numeric), \(x) coalesce(x, mean(x[idx], na.rm = TRUE))))

# modelling

mod_test <- model_matrix[-idx, ]
mod_train <- model_matrix[idx, ]

mod_train$dna_outcome <- factor(mod_train$dna_outcome, levels = c("1", "0"))
mod_test$dna_outcome <- factor(mod_test$dna_outcome, levels = c("1", "0"))

target <- "dna_outcome"
features <- names(mod_train)[!names(mod_train) %in% target]

# quick and dirty feature importance ranking

results <- data.frame(
  feature = features,
  form = NA,
  performance = NA,
  status = "Pending",
  stringsAsFactors = FALSE
)

train_samp <- mod_train[sample(seq_len(nrow(mod_train)), size = n_featimp), ]

for (i in seq_along(features)) {
  current_feature <- features[i]
  
  form_str <- paste(target, "~", current_feature)
  form <- as.formula(form_str)
  
  results$form[i] <- form_str
  
  tryCatch({
    # Run the model
    mod <- glm(form, data = train_samp, family = binomial)
    
    results$performance[i] <- summary(mod)$aic
    results$status[i] <- "Success"
  }, error = function(e) {
    results$status[i] <- paste("Error:", conditionMessage(e))
  })
}

top_feat <- results[order(results$performance), ]$feature

form <- as.formula(paste(target, "~", paste(top_feat[1:6], collapse = " + ")))

model <- ranger(
  form,
  data = mod_train,
  importance = 'permutation',
  sample.fraction = c(0.5, 0.5),
  replace = TRUE,
  probability = TRUE,
  #mtry = mtry, # Try more features per split
  min.node.size = nodes # Prevents overfitting on noise
)

# generate predictions
pred <- predict(model, data = model_matrix, type = "response")
prob_dna <- pred$predictions[, "1"]
actual <- as.numeric(as.character(model_matrix$dna_outcome))

# training performance 
  roc_data <- data.frame(actual = actual[idx], pred = prob_dna[idx])
  roc_data <- roc_data[order(-roc_data$pred), ]

  tp_cum <- cumsum(roc_data$actual == 1)
  fp_cum <- cumsum(roc_data$actual == 0)

  total_p <- sum(roc_data$actual == 1)
  total_n <- sum(roc_data$actual == 0)

  tpr <- tp_cum / total_p
  fpr <- fp_cum / total_n

  auc <- sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2);auc

















output <- dataset
