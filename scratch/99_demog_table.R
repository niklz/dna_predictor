library(gtsummary)
library(cardx)
library(dplyr)
library(forcats)
library(flextable)


model_bundle_path <- "data/processed/models/rf_final_model_500k_full.rds"
model_data_path <- "data/processed/train_engineered.rds"
model_bundle  <- readRDS(model_bundle_path)
model_data <- readRDS(model_data_path)

# ==============================================================================
# 1. BACK-TRANSFORM & LUMP EXECUTIVE FEATURES
# ==============================================================================
table1_exec_data <- model_data %>%
  mutate(
    dna_outcome = fct_recode(dna_outcome, "Attended" = "attended"),
    fct_recode(dna_outcome, "attended" = "Attended"), 
    # A. Back-transform log lead time to real days (expm1 is inverse of log1p)
    lead_time_days = expm1(lead_time_days_log),
    
    # B. Lump granular age_group factor levels into 3 broad clinical stages
    age_broad = case_when(
      age_group %in% c("<18", "18-25", "18-24") ~ "< 25 years",
      age_group %in% c("25-34", "35-44", "45-54", "55-64") ~ "25–64 years",
      age_group %in% c("65-74", "75-84", "85-94", "95+") ~ "65+ years",
      TRUE ~ "Other / Unknown"
    ),
    age_broad = factor(age_broad, levels = c("< 25 years", "25–64 years", "65+ years", "Other / Unknown")),
    
    # C. Group 10 IMD Deciles into 5 Deprivation Quintiles
    imd_clean = as.character(imd),
    imd_quintile = case_when(
      imd_clean %in% c("1", "2")  ~ "Q1 (Most Deprived)",
      imd_clean %in% c("3", "4")  ~ "Q2",
      imd_clean %in% c("5", "6")  ~ "Q3",
      imd_clean %in% c("7", "8")  ~ "Q4",
      imd_clean %in% c("9", "10") ~ "Q5 (Least Deprived)",
      TRUE                        ~ "Unknown"
    ),
    imd_quintile = factor(imd_quintile, levels = c(
      "Q1 (Most Deprived)", "Q2", "Q3", "Q4", "Q5 (Least Deprived)", "Unknown"
    )),
    
    # D. Standardise Gender categories
    gender_clean = case_when(
      gender == "Female" ~ "Female",
      gender == "Male"   ~ "Male",
      TRUE               ~ "Other / Not Stated"
    ),
    
    # E. Simplify Appointment Type
    appt_type_clean = case_when(
      grepl("NEW", appointment_type, ignore.case = TRUE) ~ "New appointment",
      TRUE                                               ~ "Follow-up appointment"
    ),
    
    # F. Collapse all 20 'a_' vulnerability/accessibility flags into a single binary flag
    has_vulnerability = if_else(
      rowSums(across(starts_with("a_"), ~ .x == 1), na.rm = TRUE) > 0, 
      "Yes", "No"
    ),
    
    # G. Recode binary flags for clean table output
    is_morning_clean = if_else(is_morning == 1, "Yes", "No"),
    has_dna_history_clean = if_else(has_dna_history == 1, "Yes", "No")
  ) %>%
  # H. Lump high-cardinality predictors (GP practice, clinic code, specialties) to Top 5 + "Other"
  mutate(across(
    any_of(c("registered_gp_practice", "clinic_code", "site_code", "local_spec_code", "national_spec_code")),
    ~ fct_lump_n(factor(.), n = 3, other_level = "Other")
  ))

# ==============================================================================
# 2. DEFINE EXECUTIVE LABELS
# ==============================================================================
table1_labels <- list(
  age_broad              ~ "Age group",
  gender_clean           ~ "Sex / Gender",
  ethnicity_group        ~ "Ethnicity group",
  imd_quintile           ~ "Deprivation quintile (IMD)",
  distance_km            ~ "Travel distance (km)",
  lead_time_days         ~ "Appointment lead time (days)",
  is_morning_clean       ~ "Morning appointment slot",
  referral_urgency       ~ "Referral urgency",
  appt_type_clean        ~ "Appointment type",
  has_vulnerability      ~ "Has accessibility / vulnerability flag",
  registered_gp_practice ~ "Registered GP practice (top 3)",
  clinic_code            ~ "Clinic code (top 3)",
  national_spec_code     ~ "Speciality (top 3)",
  has_dna_history_clean  ~ "Prior DNA history (past 12 months)"
)

# ==============================================================================
# 3. BUILD CONCISE 1-PAGE SUMMARY TABLE
# ==============================================================================
tbl_executive <- table1_exec_data %>%
  select(
    dna_outcome, age_broad, gender_clean, ethnicity_group, 
    imd_quintile, distance_km, lead_time_days, is_morning_clean, 
    referral_urgency, appt_type_clean, has_vulnerability, 
    registered_gp_practice, clinic_code, national_spec_code, has_dna_history_clean
  ) %>%
  tbl_summary(
    by = dna_outcome,
    label = table1_labels,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      lead_time_days   ~ 0 # Round lead time days to whole numbers
    ),
    missing_text = "Unknown"
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(test = list(
    all_categorical() ~ "chisq.test",
    all_continuous()  ~ "wilcox.test"
  )) %>%
  bold_labels()



# ==============================================================================
# 4. EXPORT DIRECTLY TO WORD
# ==============================================================================
tbl_executive %>%
  as_flex_table() %>%
  autofit() %>%
  save_as_docx(path = "outputs/table1_executive_summary.docx")





