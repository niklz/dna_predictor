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
  # Only remove 'Censored' patients. KEEP 'Cancelled' for the cancellation model.
  filter(final_attendance != "Censored") %>% 
  mutate(
    # METRIC 1: Wasted Capacity Outcome
    wasted_slot_outcome = ifelse(final_attendance == "DNA", 1, 0),
    
    # METRIC 2: Cancellation Behavior Outcome
    cancellation_outcome = ifelse(final_attendance == "Cancelled", 1, 0),
    
    # Standard trial setup
    trial_arm = factor(trial_arm, levels = c("Control", "Intervention", "Not in Trial")),
    
    # Per-Protocol (PP) Component Exposure Routing
    pp_exposure = case_when(
      trial_arm == "Control" ~ "1_Control", 
      trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "2_Tier2_Interactive_Only",
      trial_arm == "Intervention" & intervened_by_phone == TRUE ~ "3_Tier3_Call_Reached",
      TRUE ~ NA_character_
    )
  ) %>%
  # Convert PP exposure to factor so the regression automatically sets Control as the reference level
  mutate(
    pp_exposure = factor(pp_exposure, levels = c("1_Control", "2_Tier2_Interactive_Only", "3_Tier3_Call_Reached"))
  )


# -------------------------------------------------------------------------
# MODEL 1: Advance Cancellations (Intent-To-Treat / Arm-Level)
# -------------------------------------------------------------------------
# Evaluated on EVERYONE randomized (including DNAs and Attendances)
cancel_analysis_data <- analysis_data %>% 
  filter(trial_arm %in% c("Control", "Intervention"))

cancel_fit <- glm(cancellation_outcome ~ trial_arm + ml_baseline_risk, 
                  data = cancel_analysis_data, 
                  family = poisson(link = "log"))
cancel_model_results <- coeftest(cancel_fit, vcov = sandwich)


# -------------------------------------------------------------------------
# MODEL 2: Wasted Capacity / DNAs (Per-Protocol)
# -------------------------------------------------------------------------
# Filter out Cancellations to preserve the "At-Risk" denominator
wasted_analysis_data <- analysis_data %>% 
  filter(trial_arm %in% c("Control", "Intervention") & final_attendance != "Cancelled")

pp_wasted_fit <- glm(wasted_slot_outcome ~ pp_exposure + ml_baseline_risk, 
                     data = wasted_analysis_data, 
                     family = poisson(link = "log"))
wasted_model_results <- coeftest(pp_wasted_fit, vcov = sandwich)


# -------------------------------------------------------------------------
# Print Key Recovered Estimates
# -------------------------------------------------------------------------
print("--- MODEL 1: ADVANCE CANCELLATION (ITT) ---")
print(exp(cancel_model_results[, "Estimate"]))

print("--- MODEL 2: WASTED SLOT / DNA (PP) ---")
print(exp(wasted_model_results[, "Estimate"]))

# Validation Plotting
plot_clean_dual_patchwork(cancel_model_results, wasted_model_results, true_rr_cancel = 1.25) # Note: Make sure 1.25 matches your simulation default!


# =========================================================================
# SECTION 3: MONTE CARLO PARAMETER STABILITY LOOP
# =========================================================================

run_stability_test <- function(iterations = 50, weeks_to_simulate = 26) {
  
  map_dfr(1:iterations, function(i) {
    # 1. Simulate fresh trial
    iter_data <- simulate_clinical_trial_advanced(weeks_to_simulate = weeks_to_simulate)
    
    # 2. Base Data Prep
    iter_analysis <- iter_data %>%
      filter(final_attendance != "Censored") %>%
      mutate(
        wasted_slot_outcome = ifelse(final_attendance == "DNA", 1, 0),
        cancellation_outcome = ifelse(final_attendance == "Cancelled", 1, 0),
        trial_arm = factor(trial_arm, levels = c("Control", "Intervention")),
        pp_exposure = factor(case_when(
          trial_arm == "Control" ~ "1_Control", 
          trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "2_Tier2_Interactive_Only",
          trial_arm == "Intervention" & intervened_by_phone == TRUE ~ "3_Tier3_Call_Reached"
        ), levels = c("1_Control", "2_Tier2_Interactive_Only", "3_Tier3_Call_Reached"))
      )
    
    # 3. ITT Cancellation Model
    cancel_fit <- glm(cancellation_outcome ~ trial_arm + ml_baseline_risk, 
                      data = filter(iter_analysis, trial_arm %in% c("Control", "Intervention")), 
                      family = poisson(link = "log"))
    cancel_res <- coeftest(cancel_fit, vcov = sandwich)
    
    # 4. PP DNA Model
    dna_fit <- glm(wasted_slot_outcome ~ pp_exposure + ml_baseline_risk, 
                   data = filter(iter_analysis, trial_arm %in% c("Control", "Intervention") & final_attendance != "Cancelled"), 
                   family = poisson(link = "log"))
    dna_res <- coeftest(dna_fit, vcov = sandwich)
    
    # Extract Estimates cleanly
    tibble(
      Iteration = i,
      Cancel_ITT_Effect = exp(cancel_res["trial_armIntervention", "Estimate"]),
      DNA_Tier2_Effect  = exp(dna_res["pp_exposure2_Tier2_Interactive_Only", "Estimate"]),
      DNA_Tier3_Effect  = exp(dna_res["pp_exposure3_Tier3_Call_Reached", "Estimate"])
    )
  })
}

# Run the stability test (optional, uncomment to execute)
# stability_results <- run_stability_test(iterations = 50, weeks_to_simulate = 26)