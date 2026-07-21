library(simmer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(broom)

source("99_utils.R")

# Execute clean data population
trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = 12)
plot_trial_journeys(trial_data)

flow_volumes <- analyze_pathway_volumes(trial_data)


# Execute transition distribution plot
plot_transition_time_distributions(trial_data)

# =========================================================================
# SECTION 2: CLINICAL TRIAL EVALUATION SUITE
# =========================================================================

# Step 1: Filter Completed Manifest & Build Clean Touchpoint Variables
analysis_data <- trial_data %>%
  filter(!final_attendance %in% c("Censored/Right-Truncated", "Cancelled/Rearranged")) %>%
  mutate(
    dna_outcome = ifelse(final_attendance == "DNA", 1, 0),
    trial_arm = factor(trial_arm, levels = c("Control", "Intervention", "Not in Trial")),
    
    # GRANULAR EXPOSURE FLAGS FOR TREATMENT COMPONENT MODELING
    # Baseline group = High-Risk Control (Receives Tier 1 + Tier 2 Passive)
    received_t2_interactive = ifelse(sms_tier2_type == "Interactive Gate", 1, 0),
    received_t3_call        = ifelse(intervened_by_phone == TRUE, 1, 0)
  )

# -------------------------------------------------------------------------
# MODEL A: Intent-to-Treat (ITT) — High-Level Policy Valuation
# -------------------------------------------------------------------------
itt_model <- glm(
  dna_outcome ~ trial_arm + ml_baseline_risk,
  data = analysis_data,
  family = poisson(link = "log")
)

itt_results <- coeftest(itt_model, vcov = sandwich)

# -------------------------------------------------------------------------
# MODEL B: Per-Protocol (PP) — Granular Pathway Component Analysis
# -------------------------------------------------------------------------

# 1. Restrict and Recode with a Single Factor Structure
high_risk_trial_cohorts <- analysis_data %>% 
  filter(trial_arm != "Not in Trial") %>%
  mutate(
    # Collapse the nested logic into mutually exclusive endpoints
    pp_exposure = case_when(
      trial_arm == "Control" ~ "1_Control", 
      trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "2_Tier2_Interactive_Only",
      trial_arm == "Intervention" & intervened_by_phone == TRUE ~ "3_Tier3_Call_Reached"
    ),
    
    # Set the baseline reference group to High-Risk Controls
    pp_exposure = factor(pp_exposure, levels = c("1_Control", "2_Tier2_Interactive_Only", "3_Tier3_Call_Reached"))
  )

# 2. Fit the Modified Poisson Model
pp_component_model <- glm(
  dna_outcome ~ pp_exposure + ml_baseline_risk,
  data = high_risk_trial_cohorts,
  family = poisson(link = "log")
)

# 3. Extract Robust Standard Errors
pp_results <- coeftest(pp_component_model, vcov = sandwich)

print("--- PER-PROTOCOL (PP) EXPOSURE RESULTS ---")
print(exp(pp_results[, "Estimate"]))

# =========================================================================
# SECTION 3: RESULTS PARSING
# =========================================================================
print("--- INTENT-TO-TREAT (ITT) RESULTS ---")
print(exp(itt_results[, "Estimate"]))

print("--- PER-PROTOCOL (PP) COMPONENT EXPOSURE RESULTS ---")
print(exp(pp_results[, "Estimate"]))


library(broom)



# Example Usage: (Run this with varying target_weeks to find when power >= 0.80)
estimated_power <- calculate_simulation_power(target_weeks = 12, iterations = 50)
print(paste("Power at 12 weeks:", estimated_power))

# Execute Plot
plot_recovered_parameters(pp_results)
