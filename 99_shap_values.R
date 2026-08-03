library(fastshap)
library(shapviz)
library(tidymodels)
library(doParallel)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggnewscale)
library(legendry)

# ==============================================================================
# STEP 1: FAST SHAP COMPUTATION (ON BAKED DATA)
# ==============================================================================

# Extract model engine & recipe
raw_rf <- extract_fit_engine(model)
rec    <- extract_recipe(model)

# Pre-bake data ONCE (Fast!)
train_baked <- bake(rec, new_data = train_raw) %>% select(-any_of(target_col))

# Background and explanation samples
set.seed(123)
sample_indices <- sample(nrow(train_raw), 1000)

bg_X_baked          <- train_baked[sample(nrow(train_baked), 300), ]
explain_sample_baked <- train_baked[sample_indices, ]

# Keep matching RAW data sample for plotting labels later
explain_sample_raw   <- train_raw[sample_indices, ] %>% 
  select(-any_of(target_col)) %>%
  mutate(
    # Apply high-level grouping to raw ethnicity for reporting
    ethnicity_group = case_when(
      ethnicity %in% c("unknown", "not stated", "not known", 
                       "not collected at this time", "not set") ~ "Unknown",
      grepl("^white", ethnicity, ignore.case = TRUE) ~ "White",
      TRUE ~ "Global majority"
    )
  )

# Prediction wrapper running directly on ranger engine (No recipe overhead!)
pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions[, 1]
}

# Run fastshap in parallel
n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores)
invisible(clusterCall(cl, function(lp) .libPaths(lp), .libPaths()))
registerDoParallel(cl)

ex_global <- explain(
  object = raw_rf,
  X = bg_X_baked,
  newdata = explain_sample_baked,
  pred_wrapper = pfun,
  nsim = 50,
  adjust = TRUE,
  parallel = TRUE,
  .packages = "ranger"
)

stopCluster(cl)
registerDoSEQ()

# ==============================================================================
# STEP 1: PREPARE DATA & ASSIGN DOMAIN GROUPS
# ==============================================================================

shap_mat <- as.data.frame(ex_global)

feature_df <- as.data.frame(lapply(colnames(shap_mat), function(col) {
  if (col %in% colnames(explain_sample_raw)) {
    explain_sample_raw[[col]]
  } else {
    explain_sample_baked[[col]]
  }
}))
colnames(feature_df) <- colnames(shap_mat)

# Lump high-cardinality categoricals
high_card_cols <- c("clinic_code", "clinic_location", "site_code", "local_spec_code")

feature_df <- feature_df %>%
  mutate(across(
    any_of(high_card_cols), 
    ~ as.character(fct_lump_n(factor(.), n = 4, other_level = "Other"))
  ))

continuous_vars <- c("distance_km", "age_at_appointment", "lead_time_days_log", "appt_hour_sin", "appt_hour_cos")

# ==============================================================================
# RECODING LOOKUPS & CLEANING DICTIONARIES
# ==============================================================================

feature_lookup <- c(
  "distance_km"        = "Distance (km)",
  "age_at_appointment" = "Age at Appointment",
  "lead_time_days_log" = "Lead Time (log days)",
  "appt_hour_sin"      = "Appt Hour (Sin)",
  "appt_hour_cos"      = "Appt Hour (Cos)",
  "local_spec_code"    = "Local Specialty",
  "national_spec_code" = "National Specialty",
  "appointment_type"   = "Appointment Type",
  "gender"             = "Gender",
  "site_code"          = "Site Code",
  "appt_dow"           = "Day of Week",
  "referral_urgency"   = "Referral Urgency",
  "clinic_code"        = "Clinic Code",
  "clinic_location"    = "Clinic Location",
  "imd"                = "IMD Decile",
  "ethnicity_group"    = "Ethnicity Group",
  "appt_month_num"     = "Appointment Month",
  "lead_over_30"       = "Lead Time > 30 Days",
  "is_morning"         = "Morning Appointment",
  "has_dna_history"    = "Prior DNA History"
)

month_lookup <- setNames(month.abb, sprintf("%02d", 1:12))

# ==============================================================================
# STEP 2: PIVOT, CLEAN & ASSIGN DOMAIN CATEGORIES
# ==============================================================================

shap_long <- shap_mat %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = "feature", values_to = "shap_value")

feature_long <- feature_df %>%
  mutate(row_id = row_number()) %>%
  mutate(across(-row_id, as.character)) %>% 
  pivot_longer(-row_id, names_to = "feature", values_to = "feature_value")

full_df <- shap_long %>%
  left_join(feature_long, by = c("row_id", "feature")) %>%
  mutate(
    is_continuous = feature %in% continuous_vars,
    
    # Outer tier label: Feature Name
    feature_clean = recode(feature, !!!feature_lookup, .default = feature),
    
    # Inner tier label: Feature Value
    feature_value_clean = case_when(
      # Assign a single space to continuous variables (so they don't print "NA" or break)
      is_continuous ~ " ", 
      
      # Binary flags
      feature %in% c("is_morning", "has_dna_history", "lead_over_30") & feature_value %in% c("1", "1.0", "TRUE") ~ "Yes",
      feature %in% c("is_morning", "has_dna_history", "lead_over_30") & feature_value %in% c("0", "0.0", "FALSE") ~ "No",
      
      # Month names
      feature == "appt_month_num" & feature_value %in% names(month_lookup) ~ month_lookup[feature_value],
      feature == "appt_month_num" ~ paste0("Month ", feature_value),
      
      # Deprivation Deciles
      feature == "imd" ~ paste0("Decile ", feature_value),
      
      TRUE ~ feature_value
    ),
    
    # Domain Groups
    feature_group = case_when(
      feature %in% c("age_at_appointment", "gender", "ethnicity_group", 
                     "imd", "distance_km", "has_dna_history") ~ "Demographics & background",
      
      feature %in% c("lead_time_days_log", "lead_over_30", "appt_dow", 
                     "appt_month_num", "appt_hour_sin", "appt_hour_cos", "is_morning") ~ "Scheduling & timing",
      
      TRUE ~ "Clinical & service context"
    ),
    
    # Continuous scaling
    num_scaled = if_else(
      is_continuous,
      (suppressWarnings(as.numeric(feature_value)) - min(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE)) / 
        (max(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE) - min(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE)),
      NA_real_
    )
  )

# Set factor order for facet domain panels
full_df <- full_df %>%
  mutate(feature_group = factor(feature_group, levels = c(
    "Demographics & background",
    "Scheduling & timing",
    "Clinical & service context"
  )))

# ==============================================================================
# STEP 3: FACTOR ORDERING FOR LEGENDRY INTERACTION
# ==============================================================================

parent_order <- full_df %>%
  group_by(feature_clean) %>%
  summarise(parent_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
  arrange(parent_shap) %>%
  pull(feature_clean)

level_order <- full_df %>%
  group_by(feature_value_clean) %>%
  summarise(level_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
  arrange(level_shap) %>%
  pull(feature_value_clean)

# Setting factor levels ensures interaction(inner, outer) follows SHAP importance
full_df <- full_df %>%
  mutate(
    feature_clean       = factor(feature_clean, levels = parent_order),
    feature_value_clean = factor(feature_value_clean, levels = level_order)
  )

df_cont <- full_df %>% filter(is_continuous)
df_cat  <- full_df %>% filter(!is_continuous)

# Categorical palette
n_cat_features <- length(unique(df_cat$feature_clean))
cat_palette    <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_cat_features)

# ==============================================================================
# STEP 4: FACETED ggplot WITH LEGENDRY INTERACTION AXIS
# ==============================================================================

# ==============================================================================
# STEP 4: FACETED ggplot WITH LEGENDRY INTERACTION AXIS
# ==============================================================================

p <- ggplot(
  data = full_df,
  mapping = aes(
    x = shap_value, 
    # 1. Use the bulletproof separator in interaction
    y = interaction(feature_value_clean, feature_clean, sep = "___", drop = TRUE)
  )
) +
  geom_quasirandom(
    data = df_cont,
    aes(color = num_scaled),
    alpha = 0.6, size = 1, groupOnX = FALSE
  ) +
  scale_color_gradient(
    # name = "Continuous value",
    name = "",
    low = "#c6dbef", high = "#08306b",
    breaks = c(0, 1), labels = c("Low", "High")
  ) +
  new_scale_color() +
  geom_quasirandom(
    data = df_cat,
    aes(color = feature_clean), 
    alpha = 0.7, size = 1, groupOnX = FALSE
  ) +
  scale_color_manual(
    # name = "Categorical feature",
    name = "",
    values = cat_palette
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # facet_grid(feature_group ~ ., scales = "free_y", space = "free_y") +
  facet_wrap(vars(feature_group), nrow = 1, scales = "free") +
  
  # 2. THE FIX: Tell legendry EXACTLY how to split using 'key' instead of 'delim'
  guides(y = guide_axis_nested(key = key_range_auto(sep = "___"))) +
  
  # Remove the literal interaction(...) string from the axis title
  labs(
    y = NULL,
    x = "SHAP value (impact on predicted proability of DNA)") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8.5),
    strip.text.y = element_text(size = 10, face = "bold", angle = 270)
  )

print(p)
