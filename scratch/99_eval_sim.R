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
library(purrr)
library(patchwork)

source("99_utils.R")

# =========================================================================
# SECTION 1: SIMULATION & PROCESS MINING SUITE
# =========================================================================
set.seed(42)
trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = 26)

# Examine the operational process dynamics
plot_trial_journeys(trial_data)
process_diags <- generate_trial_process_suite(trial_data)
flow_volumes <- analyze_pathway_volumes(trial_data)
plot_transition_time_distributions(trial_data)

# =========================================================================
# SECTION 2: CLINICAL TRIAL EVALUATION SUITE
# =========================================================================

# Step 1: Base Analysis Dataset (Minimal Set)
analysis_data <- trial_data %>%
  # Only remove 'censored' patients. KEEP 'cancelled' for the cancellation model.
  filter(final_attendance != "censored") %>% 
  mutate(
    # METRIC 1: Wasted Capacity Outcome
    wasted_slot_outcome = ifelse(final_attendance == "dna", 1, 0),
    
    # METRIC 2: Cancellation Behavior Outcome
    cancellation_outcome = ifelse(final_attendance == "cancelled", 1, 0),
    
    # Standard trial setup
    trial_arm = factor(trial_arm, levels = c("control", "intervention", "not_in_trial")),
    
    # Per-Protocol (PP) Component Exposure Routing
    pp_exposure = case_when(
      trial_arm == "control" ~ "control", 
      trial_arm == "intervention" & intervened_by_phone == FALSE ~ "tier2_interactive_only",
      trial_arm == "intervention" & intervened_by_phone == TRUE ~ "tier3_call_reached",
      TRUE ~ NA_character_
    )
  ) %>%
  # Convert PP exposure to factor so the regression automatically sets control as the reference level
  mutate(
    pp_exposure = factor(pp_exposure, levels = c("control", "tier2_interactive_only", "tier3_call_reached"))
  )

# -------------------------------------------------------------------------
# MODEL 1: Advance Cancellations (Intent-To-Treat / Arm-Level)
# -------------------------------------------------------------------------
# Evaluated on EVERYONE randomized (including DNAs and Attendances)
cancel_analysis_data <- analysis_data %>% 
  filter(trial_arm %in% c("control", "intervention"))

cancel_fit <- glm(cancellation_outcome ~ trial_arm + ml_baseline_risk, 
                  data = cancel_analysis_data, 
                  family = poisson(link = "log"))
cancel_model_results <- coeftest(cancel_fit, vcov = sandwich)

# -------------------------------------------------------------------------
# MODEL 2: Wasted Capacity / DNAs (Per-Protocol)
# -------------------------------------------------------------------------
# Filter out Cancellations to preserve the "At-Risk" denominator
wasted_analysis_data <- analysis_data %>% 
  filter(trial_arm %in% c("control", "intervention") & final_attendance != "cancelled")

pp_wasted_fit <- glm(wasted_slot_outcome ~ pp_exposure + ml_baseline_risk, 
                     data = wasted_analysis_data, 
                     family = poisson(link = "log"))
wasted_model_results <- coeftest(pp_wasted_fit, vcov = sandwich)

# -------------------------------------------------------------------------
# MODEL 3: Global Clinic Capacity / DNAs (Intent-To-Treat / Arm-Level)
# -------------------------------------------------------------------------
# Evaluated on EVERYONE randomized. 
# Here, a "wasted slot" counts both DNAs AND unmanaged late disruptions, 
# or you specifically model DNA status across the whole arm.

global_itt_data <- analysis_data %>%  
  filter(trial_arm %in% c("control", "intervention"))

itt_wasted_fit <- glm(wasted_slot_outcome ~ trial_arm + ml_baseline_risk, 
                      data = global_itt_data, 
                      family = poisson(link = "log"))

itt_wasted_results <- coeftest(itt_wasted_fit, vcov = sandwich)

# -------------------------------------------------------------------------
# Print Key Recovered Estimates
# -------------------------------------------------------------------------
print("--- Model 1: advance cancellation (ITT) ---")
print(exp(cancel_model_results[, "Estimate"]))

print("--- Model 2: wasted slot / DNA (PP) ---")
print(exp(wasted_model_results[, "Estimate"]))



# =========================================================================
# SECTION 3: MONTE CARLO STABILITY
# =========================================================================
# Run the dual stability test loop
stability_results <- run_stability(iterations = 100)
plot_stability(stability_results)
