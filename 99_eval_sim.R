library(simmer)
library(bupaR)
library(processmapR)
library(processanimateR)
library(tidyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(broom)
library(dplyr)

source("99_utils.R")

# Execute clean data population
trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = 26)
# examine the process
plot_trial_journeys(trial_data)
process_diags <- generate_trial_process_suite(trial_data)
flow_volumes <- analyze_pathway_volumes(trial_data)

# Execute transition distribution plot
plot_transition_time_distributions(trial_data)

# =========================================================================
# SECTION 2: CLINICAL TRIAL EVALUATION SUITE
# =========================================================================

# Step 1: Filter Completed Manifest & Build Clean Touchpoint Variables
analysis_data <- trial_data %>%
  # FIX: Only remove 'Censored' (patients who haven't reached their appt date yet). 
  # KEEP 'Cancelled' in the dataset.
  filter(final_attendance != "Censored") %>% 
  mutate(
    # METRIC 1: The "Wasted Capacity" Outcome (Operational Perspective)
    # 1 = DNA (Total Loss), 0 = Attended OR Cancelled (Slot was utilized or recovered)
    wasted_slot_outcome = ifelse(final_attendance == "DNA", 1, 0),
    
    # METRIC 2: The "Cancellation Behavior" Outcome (Mechanism Perspective)
    # 1 = Cancelled successfully, 0 = Attended OR DNA
    cancellation_outcome = ifelse(final_attendance == "Cancelled", 1, 0),
    
    # Standard trial setup
    trial_arm = factor(trial_arm, levels = c("Control", "Intervention", "Not in Trial")),
    
    received_t2_interactive = ifelse(sms_tier2_type == "Interactive Gate", 1, 0),
    received_t3_call        = ifelse(intervened_by_phone == TRUE, 1, 0)
  )

# -------------------------------------------------------------------------
# MODEL A: Intent-to-Treat (ITT) 
# -------------------------------------------------------------------------
itt_model <- glm(
  wasted_slot_outcome ~ trial_arm + offset(log(ml_baseline_risk)),
  data = analysis_data,
  family = poisson(link = "log")
)

itt_results <- coeftest(itt_model, vcov = sandwich)

# =========================================================================
# MODEL B: Per-Protocol (PP) Component Analysis (Behaviorally Confounded)
# =========================================================================

high_risk_trial_cohorts <- analysis_data %>% 
  filter(trial_arm != "Not in Trial") %>%
  mutate(
    pp_exposure = case_when(
      trial_arm == "Control" ~ "1_Control", 
      
      # Intervention text-only group: Must not have been called
      trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "2_Tier2_Interactive_Only",
      
      # Intervention phone call group: Reached successfully by phone
      trial_arm == "Intervention" & intervened_by_phone == TRUE ~ "3_Tier3_Call_Reached",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pp_exposure)) %>%
  mutate(
    pp_exposure = factor(pp_exposure, levels = c("1_Control", "2_Tier2_Interactive_Only", "3_Tier3_Call_Reached"))
  )

# -------------------------------------------------------------------------
# 1. METRIC 1: Wasted Slot Outcome (Operational Capacity Loss)
# -------------------------------------------------------------------------
# Offset remains mathematically required because ML baseline risk 
# scales directly with overall DNA / wasted slot probability.
pp_wasted_model <- glm(
  wasted_slot_outcome ~ pp_exposure + offset(log(ml_baseline_risk)),
  data = high_risk_trial_cohorts,
  family = poisson(link = "log")
)

pp_wasted_results <- coeftest(pp_wasted_model, vcov = sandwich)


# -------------------------------------------------------------------------
# 2. METRIC 2: Cancellation Behavior (Mechanism of Action)
# -------------------------------------------------------------------------
# CRITICAL CHANGE: Because your simulation ties cancellation probability 
# to a quadratic curve of ml_baseline_risk, you MUST include 
# ml_baseline_risk as a covariate here to isolate the true intervention effect.
# (Note: No offset here, because cancellations are not a direct 1:1 scale of ML DNA risk).
pp_cancel_model <- glm(
  cancellation_outcome ~ pp_exposure + ml_baseline_risk, 
  data = high_risk_trial_cohorts,
  family = poisson(link = "log")
)

pp_cancel_results <- coeftest(pp_cancel_model, vcov = sandwich)


# Print Results
print("--- PP EXPOSURE: WASTED SLOT (Operational Loss) ---")
print(exp(pp_wasted_results[, "Estimate"]))

print("--- PP EXPOSURE: SUCCESSFUL ADVANCE CANCELLATION ---")
print(exp(pp_cancel_results[, "Estimate"]))

print("--- PER-PROTOCOL (PP) EXPOSURE RESULTS ---")
print(exp(pp_results[, "Estimate"]))

plot_recovered_parameters(pp_wasted_results)

# stability test
stability_results <- run_stability_test(iterations = 50, weeks_to_simulate = 26)





