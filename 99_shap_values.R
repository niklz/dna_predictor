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
  predict(object, data = newdata)$predictions[, 2]
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
# STEP 2: PIVOT & ASSIGN DOMAIN CATEGORIES
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
    
    # Assign Domain Groups
    feature_group = case_when(
      feature %in% c("age_at_appointment", "gender", "ethnicity_group", 
                     "imd", "distance_km", "has_dna_history") ~ "Demographics & Background",
      
      feature %in% c("lead_time_days_log", "lead_over_30", "appt_dow", 
                     "appt_month_num", "appt_hour_sin", "appt_hour_cos", "is_morning") ~ "Scheduling & Timing",
      
      TRUE ~ "Clinical & Service Context"
    ),
    
    # Continuous scaling
    num_scaled = if_else(
      is_continuous,
      (suppressWarnings(as.numeric(feature_value)) - min(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE)) / 
        (max(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE) - min(suppressWarnings(as.numeric(feature_value)), na.rm = TRUE)),
      NA_real_
    ),
    
    plot_feature = if_else(
      is_continuous, 
      feature, 
      paste0(feature, " = ", feature_value)
    )
  )

# Set logical factor order for panels
full_df <- full_df %>%
  mutate(feature_group = factor(feature_group, levels = c(
    "Demographics & Background",
    "Scheduling & Timing",
    "Clinical & Service Context"
  )))

# ==============================================================================
# STEP 3: INTERLEAVED FACTOR ORDERING WITHIN DOMAINS
# ==============================================================================

parent_importance <- full_df %>%
  group_by(feature) %>%
  summarise(parent_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop")

level_importance <- full_df %>%
  group_by(feature, plot_feature) %>%
  summarise(level_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop")

feature_order <- level_importance %>%
  left_join(parent_importance, by = "feature") %>%
  arrange(parent_shap, level_shap) %>%
  pull(plot_feature) %>%
  unique()

full_df <- full_df %>%
  mutate(plot_feature = factor(plot_feature, levels = feature_order))

df_cont <- full_df %>% filter(is_continuous)
df_cat  <- full_df %>% filter(!is_continuous)

# Dynamic palette for categoricals
n_cat_features <- length(unique(df_cat$feature))
cat_palette    <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_cat_features)

# ==============================================================================
# STEP 4: FACETED ggplot
# ==============================================================================

p <- ggplot(mapping = aes(x = shap_value, y = plot_feature)) +
  
  # Layer 1: Continuous Features
  geom_quasirandom(
    data = df_cont,
    aes(color = num_scaled),
    alpha = 0.6, size = 1, groupOnX = FALSE
  ) +
  scale_color_gradient(
    name = "Continuous Value",
    low = "#c6dbef",  
    high = "#08306b", 
    breaks = c(0, 1),
    labels = c("Low", "High")
  ) +
  
  new_scale_color() +
  
  # Layer 2: Categoricals
  geom_quasirandom(
    data = df_cat,
    aes(color = feature), 
    alpha = 0.7, size = 1, groupOnX = FALSE
  ) +
  scale_color_manual(
    name = "Categorical Feature",
    values = cat_palette
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # FACETING BY DOMAIN WITH PROPORTIONAL HEIGHTS
  # facet_grid(feature_group ~ ., scales = "free_y", space = "free_y") +
  facet_wrap(~ feature_group, nrow = 1, scales = "free_y") +
  
  labs(
    title = "Global SHAP Summary by Domain",
    subtitle = "Interleaved feature effects grouped into clinical, scheduling, and demographic categories",
    x = "SHAP Value (Impact on Prediction)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8.5),
    strip.text.y = element_text(size = 10, face = "bold", angle = 270),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing = unit(0.8, "lines")
  )

print(p)
