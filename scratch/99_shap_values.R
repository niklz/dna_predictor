library(fastshap)
library(shapviz)
library(doParallel)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggnewscale)
library(legendry)

# -------------------------------------------------------------------------
# 1. Setup & Configuration
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()

model_bundle <- readRDS("data/processed/models/rf_final_model_500k_lean.rds") 
model <- model_bundle$model
train_raw <- readRDS("data/processed/archive/train_engineered.rds")

# ==============================================================================
# STEP 1: FAST SHAP COMPUTATION (ON BAKED DATA)
# ==============================================================================


raw_rf <- extract_fit_engine(model)
rec    <- extract_recipe(model)


train_baked <- bake(rec, new_data = train_raw) %>% 
  select(
    -any_of(conf$target_col), 
    -any_of(c("dim_patient_id", "patient_profile"))
  )

# ------------------------------------------------------------------------------
# THE SHAP CAVEAT: BACKGROUND VS. EXPLANATION DATA
# ------------------------------------------------------------------------------
# 1. Background Data: MUST be a simple random sample. 
#    This sets the "base value" (average prediction). If you skew this, 
#    you artificially change what the model considers "normal".

set.seed(123)
bg_X_baked <- train_baked[sample(nrow(train_baked), 1000), ]

# 2. Explanation Data: Use Inverse Frequency Weighting to oversample rare types.
#    Pick the columns that define your sub-populations (e.g., ethnicity & urgency)
stratify_cols <- c("ethnicity_group", "referral_urgency")

# Calculate weights: 1 / (number of patients in that specific group)
weighting_df <- train_raw %>%
  select(all_of(stratify_cols)) %>%
  add_count(across(all_of(stratify_cols)), name = "group_n") %>%
  mutate(sample_weight = 1 / group_n)

# Sample indices using probabilities based on our calculated weights
set.seed(123)
sample_indices <- sample(
  seq_len(nrow(train_raw)), 
  size = 1000, 
  prob = weighting_df$sample_weight, 
  replace = FALSE 
)

explain_sample_baked <- train_baked[sample_indices, ]

# Keep matching RAW data sample for plotting labels later
explain_sample_raw   <- train_raw[sample_indices, ] %>% 
  select(-any_of(conf$target_col)) %>%
  mutate(
    ethnicity_group = case_when(
      ethnicity_group %in% c("unknown", "not stated", "not known", 
                             "not collected at this time", "not set") ~ "Unknown",
      grepl("^white", ethnicity_group, ignore.case = TRUE) ~ "White",
      TRUE ~ "Global majority"
    )
  )


pfun <- function(object, newdata) { 
  predict(object, data = newdata, num.threads = 24)$predictions[, 1] 
}

# 4. Sequential fastshap call (Bypasses Windows socket cloning)

ex_global <- fastshap::explain(
  object       = raw_rf,
  X            = bg_X_baked,
  newdata      = explain_sample_baked,
  pred_wrapper = pfun,
  nsim         = 50,         
  adjust       = TRUE,
  parallel     = FALSE,  # FALSE: Tells fastshap not to spawn messy background R sessions
  .packages    = "ranger"
)

gc()
message("-> [COMPLETE] SHAP values generated")


shap_plot <- local({
  # ==============================================================================
  # RECODING LOOKUPS & CLEANING DICTIONARIES [cite: 636]
  # ==============================================================================
  feature_lookup <- c(
    "distance_km"            = "Distance (km)",
    "age_at_appointment"     = "Age at Appointment",
    "age_group"              = "Age Group",
    "lead_time_days_log"     = "Lead Time (log days)",
    "appt_hour_sin"          = "Appt Hour (Sin)",
    "appt_hour_cos"          = "Appt Hour (Cos)",
    "local_spec_code"        = "Local Specialty",
    "national_spec_code"     = "National Specialty",      # Recodes national specialty
    "appointment_type"       = "Appointment Type",
    "gender"                 = "Gender",
    "site_code"              = "Site Code",
    "appt_dow"               = "Day of Week",
    "referral_urgency"       = "Referral Urgency",
    "clinic_code"            = "Clinic Code",
    "clinic_location"        = "Clinic Location",
    "imd"                    = "IMD Decile",
    "ethnicity_group"        = "Ethnicity Group",
    "appt_month_num"         = "Appointment Month",
    "lead_over_30"           = "Lead Time > 30 Days",
    "is_morning"             = "Morning Appointment",
    "has_dna_history"        = "Prior DNA History",
    "registered_gp_practice" = "GP Practice"              # Recodes GP practice
  )
  
  month_lookup <- setNames(month.abb, sprintf("%02d", 1:12))
  
  # ==============================================================================
  # STEP 1: PREPARE DATA & DOMAIN GROUPS [cite: 637]
  # ==============================================================================
  shap_mat <- as.data.frame(ex_global)
  
  # Explicit list of continuous variables to handle separately [cite: 637]
  continuous_vars <- c("distance_km", "age_at_appointment", "lead_time_days_log", "appt_hour_sin", "appt_hour_cos")
  
  # FIXED: Added registered_gp_practice and national_spec_code to back-transformation [cite: 637]
  vars_to_back_transform <- c(
    "clinic_location", "site_code", "local_spec_code", "appointment_type", 
    "registered_gp_practice", "national_spec_code"
  )
  
  feature_df <- as.data.frame(lapply(colnames(shap_mat), function(col) { 
    if (col %in% vars_to_back_transform && col %in% colnames(explain_sample_raw)) { 
      explain_sample_raw[[col]] 
    } else { 
      explain_sample_baked[[col]] 
    } 
  })) 
  colnames(feature_df) <- colnames(shap_mat)
  
  # FIXED: Added registered_gp_practice and national_spec_code to lumping logic [cite: 638]
  high_card_cols <- c(
    "clinic_code", "clinic_location", "site_code", "local_spec_code", 
    "registered_gp_practice", "national_spec_code"
  ) 
  feature_df <- feature_df %>% mutate(across(
    any_of(high_card_cols),
    ~ as.character(forcats::fct_lump_n(factor(.), n = 5, other_level = "Other"))
  ))
  
  # ==============================================================================
  # STEP 2: PIVOT, CLEAN & ASSIGN DOMAIN CATEGORIES [cite: 639]
  # ==============================================================================
  shap_long <- shap_mat %>% 
    mutate(row_id = row_number()) %>% 
    pivot_longer(-row_id, names_to = "feature", values_to = "shap_value")
  
  feature_long <- feature_df %>% 
    mutate(row_id = row_number()) %>% 
    mutate(across(-row_id, as.character)) %>% 
    pivot_longer(-row_id, names_to = "feature", values_to = "feature_value")
  
  full_df <- left_join(shap_long, feature_long, by = c("row_id", "feature")) %>% 
    mutate(
      feature_clean = recode(feature, !!!feature_lookup, .default = feature),
      is_continuous = feature %in% continuous_vars,
      
      # Scale continuous features for gradient colors [cite: 640]
      num_scaled = if_else(
        is_continuous,
        (suppressWarnings(as.numeric(feature_value)) - min(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE)) / 
          (max(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE) - min(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE)),
        NA_real_
      )
    )
  
  # Clean up categorical values [cite: 641]
  full_df <- full_df %>% 
    mutate(
      feature_value_clean = case_when(
        # Collapse continuous variables to prevent discrete coordinate breaks [cite: 641]
        is_continuous ~ " ",
        
        # Binary flags [cite: 641]
        feature %in% c("is_morning", "has_dna_history", "lead_over_30") & feature_value %in% c("1", "1.0", "TRUE") ~ "Yes",
        feature %in% c("is_morning", "has_dna_history", "lead_over_30") & feature_value %in% c("0", "0.0", "FALSE") ~ "No",
        
        # Month names [cite: 641]
        feature == "appt_month_num" & feature_value %in% names(month_lookup) ~ month_lookup[feature_value],
        feature == "appt_month_num" ~ paste0("Month ", feature_value),
        
        # Deprivation Deciles [cite: 641]
        feature == "imd" ~ paste0("Decile ", feature_value),
        
        # TWEAK 1: Standardise "other"/"Other" categories across all features [cite: 641]
        tolower(feature_value) %in% c("other", "Decile other") ~ "Other",
        
        TRUE ~ as.character(feature_value)
      )
    )
  
  # FIXED: mapped has_dna_history to "Clinical & service context" (by omitting it from demographics) [cite: 642]
  full_df <- full_df %>% 
    mutate(
      feature_group = case_when(
        feature %in% c("has_dna_history", "age_group", "age_at_appointment", "gender", "ethnicity_group", "imd", "distance_km") ~ "Demographics & background",
        feature %in% c("appt_dow", "appt_month_num", "lead_time_days_log", "lead_over_30", "is_morning", "appt_hour_sin", "appt_hour_cos") ~ "Scheduling & timing",
        TRUE ~ "Clinical & service context"
      )
    ) %>% 
    mutate(
      feature_group = factor(feature_group, levels = c(
        "Demographics & background", "Scheduling & timing", "Clinical & service context"
      ))
    )
  
  # ==============================================================================
  # STEP 3: DYNAMIC FACTOR LEVEL EXTRACTION & RELEVELING [cite: 643]
  # ==============================================================================
  if (!"age_group" %in% names(feature_lookup)) { 
    feature_lookup <- c(feature_lookup, "age_group" = "Age Group") 
  }
  
  # Determine outer y-axis sorting (absolute feature impact) [cite: 643]
  parent_order <- full_df %>% 
    group_by(feature_clean) %>% 
    summarise(parent_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>% 
    arrange(parent_shap) %>% 
    pull(feature_clean)
  
  # Automatically extract pre-existing factor levels from source datasets [cite: 644]
  factor_cols <- c() 
  if (exists("train_baked")) factor_cols <- unique(c(factor_cols, colnames(train_baked)[sapply(train_baked, is.factor)])) 
  if (exists("train_raw"))   factor_cols <- unique(c(factor_cols, colnames(train_raw)[sapply(train_raw, is.factor)]))
  
  natural_features_levels <- list() 
  for (col in factor_cols) { 
    clean_parent_name <- recode(col, !!!feature_lookup, .default = col)
    
    raw_levels <- if (exists("train_baked") && col %in% colnames(train_baked) && is.factor(train_baked[[col]])) {
      levels(train_baked[[col]])
    } else if (exists("train_raw") && col %in% colnames(train_raw) && is.factor(train_raw[[col]])) {
      levels(train_raw[[col]])
    } else {
      NULL
    }
    
    if (!is.null(raw_levels)) {
      clean_levels <- sapply(raw_levels, function(val) {
        case_when(
          col %in% c("is_morning", "has_dna_history", "lead_over_30") && val %in% c("1", "1.0", "TRUE") ~ "Yes",
          col %in% c("is_morning", "has_dna_history", "lead_over_30") && val %in% c("0", "0.0", "FALSE") ~ "No",
          col == "appt_month_num" && val %in% names(month_lookup) ~ month_lookup[val],
          col == "appt_month_num" ~ paste0("Month ", val),
          col == "imd" ~ paste0("Decile ", val),
          
          # TWEAK 2: Standardise natural factor levels (forces "other" -> "Other") [cite: 644]
          tolower(val) %in% c("other", "Decile other") ~ "Other",
          
          TRUE ~ as.character(val)
        )
      })
      natural_features_levels[[clean_parent_name]] <- as.character(clean_levels)
    }
  }
  
  # Sort standard nominal features by SHAP importance [cite: 645]
  natural_parent_names <- names(natural_features_levels)
  shap_sorted_values <- full_df %>% 
    filter(!feature_clean %in% natural_parent_names) %>% 
    group_by(feature_value_clean) %>% 
    summarise(level_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>% 
    arrange(level_shap) %>% 
    pull(feature_value_clean)
  
  shap_sorted_values <- shap_sorted_values[!shap_sorted_values %in% c("Other", "other")]
  
  unlisted_natural_levels <- unlist(natural_features_levels, use.names = FALSE) 
  all_present_values      <- unique(full_df$feature_value_clean) 
  unrepresented_values    <- setdiff(all_present_values, c(shap_sorted_values, unlisted_natural_levels)) 
  
  final_level_order <- unique(c(shap_sorted_values, unlisted_natural_levels, unrepresented_values)) 
  
  # Apply factor levels safely 
  full_df <- full_df %>% 
    mutate( 
      feature_clean       = factor(feature_clean, levels = parent_order), 
      feature_value_clean = factor(feature_value_clean, levels = final_level_order) 
    ) %>% 
    mutate( 
      # Move "Other"/"other" categories to the absolute bottom of each feature list 
      feature_value_clean = forcats::fct_relevel(feature_value_clean, "Other", "other", after = 0) 
    )
  
  # Split datasets for continuous vs categorical geom mapping 
  df_cont <- full_df %>% filter(is_continuous) 
  df_cat  <- full_df %>% filter(!is_continuous) 
  
  # ==============================================================================
  # STEP 4: FACETED GGPLOT WITH LEGENDRY INTERACTION AXIS 
  # ==============================================================================
  n_cat_features <- length(unique(df_cat$feature_clean))
  cat_palette    <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_cat_features) 
  
  p <- ggplot( 
    data = full_df, 
    mapping = aes( 
      x = shap_value, 
      y = interaction(feature_value_clean, feature_clean, sep = "___", drop = TRUE) 
    ) 
  ) + 
    # Continuous variables (gradient colors) [cite: 647] 
    geom_quasirandom( 
      data = df_cont, 
      aes(color = num_scaled), 
      alpha = 0.6, 
      size = 1, 
      groupOnX = FALSE 
    ) + 
    scale_color_gradient( 
      name = "", 
      low = "#c6dbef", 
      high = "#08306b", 
      breaks = c(0, 1), 
      labels = c("Low", "High") 
    ) + 
    new_scale_color() +
    # Categorical variables [cite: 648]
    geom_quasirandom(
      data = df_cat,
      aes(color = feature_clean),
      alpha = 0.7,
      size = 1,
      groupOnX = FALSE
    ) +
    scale_color_manual(
      name = "",
      values = cat_palette
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Panel layout split by Domain Group [cite: 648]
    facet_wrap(vars(feature_group), nrow = 1, scales = "free") +
    
    # Nested interaction separator [cite: 648]
    guides(y = guide_axis_nested(key = key_range_auto(sep = "___"))) +
    
    labs(
      y = NULL,
      x = "SHAP value (impact on predicted probability of DNA)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8.5),
      strip.text = element_text(size = 10, face = "bold")
    )


  p
});print(shap_plot)

