# --- Load Data ---
dataset <- read.csv("data/dna_combined.csv")


library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(stringr)
library(forcats)
library(ranger)


# --- Hyperparams ---
n_lump <- 10
n_featimp <- 5E3
nodes <- 250
mtry <- 10

# --- Initial Pre-processing ---
dataset$dna_outcome <- factor(ifelse(dataset$AttendedStatusCode == "3", "DNA", "attended"))
dataset$imd <- as.character(dataset$Index_Multiple_Deprivation_Decile)
dataset$imd <- ifelse(is.na(dataset$imd), "unknown", dataset$imd)

target_col <- "dna_outcome"
vars <- c(
  "Local_Spec_code", "National_Spec_Code", "AppointmentType", "distance_km",
  "NFA_Ind", "Age.Group", "AgeAtAppointment", "Ethnicity", "a_LD",
  "a_Autism", "a_InterpreterReq_BSL", "a_InterpreterReq_lang", "a_Balance",
  "a_CognitiveImpairment", "a_MobilityRestriction", "a_HearVisImpaired",
  "a_Dementia", "a_Depression", "a_DownsSyndrome", "a_LongStandingCondition",
  "a_Makaton", "a_MildCognitiveImpairment", "a_MemoryImpairment", "a_MoodDisorder",
  "a_OtherDisability", "a_Psychosis", "a_SevereAnxiety", "a_WheelchairUser",
  "Gender", "RegisteredGPPractice", "SITE_CODE", "prev_dna_LY",
  "Appt_Hour", "Appt_DOW", "Appt_Month", "Appt_Wknd_Ind",
  "Referral_Urgency", "Lead_Time_days", "CLINIC_CODE", "CLINIC_LOCATION", "imd"
)

# Feature engineering that doesn't depend on target/global stats
process_features <- function(df) {
  df %>%
    select(any_of(vars)) %>%
    rename_with(str_to_lower) %>%
    mutate(across(where(is.character), factor)) %>%
    mutate(
      neurodiv = if_any(starts_with("a_"), ~ .x == 1),
      neurodiv_count = rowSums(select(., starts_with("a_"))),
      appt_month = factor(substring(appt_month, 1, 10) %>% as.Date(format = "%d/%m/%Y") %>% format("%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      is_school_run = ifelse(appt_hour %in% c(8, 9, 15, 16), 1, 0),
      has_dna_history = prev_dna_ly > 0,
      is_local = distance_km < 5, 
      is_very_far = distance_km > 20
    )
}

# --- 1. Split Data into Training and Test ---
idx_train <- which(dataset$test_train == "Training")

train_raw <- dataset[idx_train, ]
test_raw  <- dataset[-idx_train, ]

train_feat <- process_features(train_raw)
test_feat  <- process_features(test_raw)

# --- 2. Define Encoders (Train Only) ---
high_card_cols <- train_feat %>%
  select(where(is.factor)) %>%
  summarise(across(everything(), n_distinct)) %>%
  pivot_longer(cols = everything()) %>%
  filter(value > 20) %>%
  pull(name)

train_lmm_encoders <- function(data, target_col, high_card_cols) {
  map(high_card_cols, function(col) {
    form <- as.formula(paste0(target_col, " ~ (1 | ", col, ")"))
    fit <- glmer(form, data = data, family = binomial)
    effects <- ranef(fit)[[1]]
    colnames(effects) <- "rel_adj"
    effects[[col]] <- rownames(effects)
    list(mapping = effects, global_intercept = fixef(fit)[[1]], col_name = col)
  }) %>% set_names(high_card_cols)
}

# Apply Encoders helper
apply_lmm <- function(x, col_name, encoder_list) {
  enc <- encoder_list[[col_name]]
  res <- enc$mapping$rel_adj[match(x, enc$mapping[[col_name]])]
  enc$global_intercept + ifelse(is.na(res), 0, res)
}

# Generate encoders using only training data
my_encoders <- train_lmm_encoders(bind_cols(target = train_raw$dna_outcome, train_feat), "target", high_card_cols)

# Calculate numeric means using only training data
train_means <- train_feat %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

# --- 3. Unified Transformation Pipeline ---
apply_pipeline <- function(feat_df, target_vec, encoders, means, high_card) {
  processed <- feat_df %>%
    # Encode high cardinality
    mutate(across(all_of(high_card), ~ apply_lmm(.x, cur_column(), encoders), .names = "{.col}_encoded")) %>%
    select(!any_of(high_card)) %>%
    # Impute factors
    mutate(across(where(is.factor), \(x) fct_na_value_to_level(x, level = "missing")))
  
  # Impute numerics using pre-calculated training means
  for(col in names(means)) {
    if(col %in% names(processed)) {
      processed[[col]][is.na(processed[[col]])] <- means[[col]]
    }
  }
  
  processed$dna_outcome <- target_vec
  return(processed)
}

mod_train <- apply_pipeline(train_feat, train_raw$dna_outcome, my_encoders, train_means, high_card_cols)
mod_test  <- apply_pipeline(test_feat, test_raw$dna_outcome, my_encoders, train_means, high_card_cols)

# --- 4. Feature Selection ---
features_list <- names(mod_train)[names(mod_train) != "dna_outcome"]
train_samp <- mod_train[sample(seq_len(nrow(mod_train)), size = min(nrow(mod_train), n_featimp)), ]

results <- map_df(features_list, function(f) {
  tryCatch({
    mod <- glm(as.formula(paste("dna_outcome ~", f)), data = train_samp, family = binomial)
    data.frame(feature = f, performance = summary(mod)$aic, status = "Success")
  }, error = function(e) data.frame(feature = f, performance = NA, status = "Error"))
})

top_feat <- results %>% filter(status == "Success") %>% arrange(performance) %>% pull(feature)
final_form <- as.formula(paste("dna_outcome ~", paste(top_feat, collapse = " + ")))

# --- 5. Modelling ---
model <- ranger(
  final_form,
  data = mod_train,
  importance = 'permutation',
  sample.fraction = c(0.5, 0.5),
  mtry = mtry,
  replace = TRUE,
  probability = TRUE,
  min.node.size = nodes
)

# --- 6. Performance Evaluation ---
calc_auc <- function(actual, predicted_prob) {
  df <- data.frame(actual = actual, pred = predicted_prob) %>% arrange(-pred)
  tp_cum <- cumsum(df$actual == "DNA")
  fp_cum <- cumsum(df$actual == "attended")
  tpr <- tp_cum / sum(df$actual == "DNA")
  fpr <- fp_cum / sum(df$actual == "attended")
  sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2)
}

# OOB Performance
train_pred_oob <- model$predictions[, "DNA"]
auc_train_oob <- calc_auc(mod_train$dna_outcome, train_pred_oob)


output <- dataset