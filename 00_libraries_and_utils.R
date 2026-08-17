library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(stringr)
library(forcats)
library(ranger)
library(ggplot2)
library(vip)
library(tidymodels)
library(probably)
library(embed) # For step_lmer (High-cardinality encoding)
library(themis)
library(tictoc)
library(future)
library(doFuture)
library(GGally)
library(config)



# -------------------------------------------------------------------------
# CUSTOM HELPER FUNCTIONS
# -------------------------------------------------------------------------
collect_ethnicity <- function(x) {
  dplyr::case_when(
    x %in% c("Not Stated", "Not Known", "Not Collected At This Time", "Not Set") ~ "Unknown",
    x %in% c("White - British", "White - Any other White background", "White - Irish", "Any other White background") ~ "White",
    x %in% c("Asian or Asian British - Indian", "Asian - Indian or British Indian", 
             "Asian or Asian British - Pakistani", "Asian - Pakistani or British Pakistani",
             "Asian or Asian British - Bangladeshi", "Asian - Bangladeshi or British Bangladeshi",
             "Asian or Asian British - Any other Asian background", "Asian - Any other Asian background") ~ "Asian",
    x %in% c("Black or Black British - Caribbean", "Black - Caribbean or Black British Caribbean",
             "Black Caribbean", "Black or Black British - African", 
             "Black  - African or British African", "Black or Black British - Any other Black background",
             "Black - Any other Black background") ~ "Black",
    x %in% c("Mixed - White and Black African", "Mixed - Any other background", 
             "Mixed - White and Black Caribbean", "Mixed - White and Asian", "Any other Mixed Background") ~ "Mixed",
    x %in% c("Other Ethnic Group - Chinese", "Chinese", "Any Other Ethnic Group", "Any other Ethnic Group") ~ "Other",
    TRUE ~ "Unknown"
  )
}

aggregate_ethnicity_high_level <- function(x) {
  x_clean <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_clean %in% c("unknown", "not stated", "not known", "not collected at this time", "not set") ~ "Unknown",
    grepl("^white", x_clean) ~ "White",
    grepl("^asian", x_clean) ~ "Asian",
    grepl("^black", x_clean) ~ "Black",
    grepl("^mixed", x_clean) ~ "Mixed",
    grepl("chinese", x_clean) | grepl("other", x_clean) ~ "Chinese",
    TRUE ~ "Unknown"
  )
}

apply_custom_feature_engineering <- function(data) {
  data %>%
    mutate(
      ethnicity_clean = collect_ethnicity(ethnicity),
      ethnicity_group = aggregate_ethnicity_high_level(ethnicity_clean),
      # Set the baseline reference level to the most frequent category
      ethnicity_group = relevel(factor(ethnicity_group), ref = "White"),
      appt_date = as.Date(substring(appt_month, 1, 10), format = "%d/%m/%Y"),
      appt_dow = factor(weekdays(appt_date)),
      appt_month_num = as.factor(format(appt_date, "%m")),
      lead_over_30 = ifelse(lead_time_days > 30, 1, 0),
      lead_time_days_log = log1p(pmax(0, lead_time_days)),
      is_morning = ifelse(appt_hour < 12, 1, 0),
      appt_hour_sin = sin(2 * pi * appt_hour / 24),
      appt_hour_cos = cos(2 * pi * appt_hour / 24),
      has_dna_history = ifelse(prev_dna_ly > 0, 1, 0)
    ) %>%
    # Replaces step_rm()
    select(
      -appt_hour, -lead_time_days, -appt_date, -appt_month, 
      -prev_dna_ly, -ethnicity, -ethnicity_clean, -age_group
    )
}
