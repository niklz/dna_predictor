library(simmer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sandwich)
library(lmtest)

# =========================================================================
# SECTION 1: DATA GENERATION ENGINE
# =========================================================================
simulate_clinical_trial_advanced <- function(
    weeks_to_simulate = 8,
    total_patients_per_week = 1000,
    ml_high_risk_prop = 0.25,
    trial_allocation_ratio = 0.50,
    
    # Coordinator Constraints
    coordinator_capacity = 1,
    max_calls_per_week = 100,
    response_prob = 0.30,
    prob_wish_to_cancel = 0.20,
    
    # INCREMENTAL EFFECT SIZES (Risk Ratios)
    rr_tier1          = 0.95,  # 5% relative reduction from universal text
    rr_tier2_passive  = 0.90,  # 10% relative reduction from passive text
    rr_tier2_interact = 0.85,  # 15% relative reduction from interactive gate
    rr_tier3_call     = 0.40   # 60% relative reduction if reached by phone
) {
  
  mins_per_day <- 1440
  mins_per_week <- 7 * mins_per_day
  total_sim_time <- weeks_to_simulate * mins_per_week
  work_mins_per_week <- 5 * 8 * 60
  mean_handling_time <- work_mins_per_week / max_calls_per_week
  
  # Step 1: Cohort Setup & ML Baseline Risk Assignment
  num_weeks <- weeks_to_simulate
  weekly_total_sizes <- rpois(num_weeks, lambda = total_patients_per_week)
  
  trial_manifest <- data.frame(week = rep(1:num_weeks, times = weekly_total_sizes)) %>%
    mutate(
      patient_id = paste0("PT_", row_number()),
      risk_profile = ifelse(runif(n()) < ml_high_risk_prop, "High-Risk", "Low-Risk"),
      
      ml_baseline_risk = ifelse(risk_profile == "Low-Risk", 
                                rbeta(n(), shape1 = 2, shape2 = 50), 
                                rbeta(n(), shape1 = 5, shape2 = 20)),
      
      trial_arm = case_when(
        risk_profile == "Low-Risk" ~ "Not in Trial",
        risk_profile == "High-Risk" & runif(n()) < trial_allocation_ratio ~ "Control",
        TRUE ~ "Intervention"
      )
    )
  
  sim_batch_sizes <- trial_manifest %>% filter(trial_arm == "Intervention") %>% count(week) %>% pull(n)
  weekly_timestamps <- seq(0, by = mins_per_week, length.out = num_weeks)
  batch_arrival_times <- rep(weekly_timestamps, times = sim_batch_sizes)
  
  # Step 2: DES Simulation Process
  env <- simmer("trial_pathway")
  patient_path <- trajectory("High_Risk_Intervention_Journey") %>%
    set_attribute("arrival_time", function() { now(env) }) %>%
    set_attribute("days_until_appt", function() { runif(1, min = 4, max = 14) }) %>%
    set_attribute("appt_time", function() { now(env) + (get_attribute(env, "days_until_appt") * mins_per_day) }) %>%
    set_attribute("tier1_text_time", function() { get_attribute(env, "appt_time") - (7 * mins_per_day) }) %>%
    set_attribute("responds_to_text", function() { rbinom(1, 1, response_prob) }) %>%
    
    timeout(function() { max(0, (get_attribute(env, "days_until_appt") - 3) * mins_per_day) }) %>%
    set_attribute("tier2_text_time", function() { now(env) }) %>%
    
    renege_in(
      t = function() { max(0, get_attribute(env, "appt_time") - now(env)) },
      out = trajectory() %>% set_attribute("outcome_status", 0) 
    ) %>%
    seize("coordinator", 1) %>%
    renege_abort() %>% 
    
    set_attribute("call_time", function() { now(env) }) %>%
    timeout(function() { max(5, rnorm(1, mean = mean_handling_time, sd = 5)) }) %>%
    set_attribute("wants_to_cancel", function() { rbinom(1, 1, prob_wish_to_cancel) }) %>%
    release("coordinator", 1) %>%
    set_attribute("outcome_status", 1)
  
  working_hours <- schedule(timetable = c(0, 9*60, 17*60), values = c(0, coordinator_capacity, 0), period = mins_per_day)
  
  env %>%
    add_resource("coordinator", capacity = working_hours) %>%
    add_generator("Intervention_Patient_", patient_path, at(batch_arrival_times), mon = 2) %>%
    run(until = total_sim_time)
  
  # Step 3: Extract and Safely Flatten Attributes
  # 1. Extract the simmer identity index. 
  # simmer is 0-indexed, so we add 1 to match R's 1-based row_number()
  sim_arrivals <- get_mon_arrivals(env) %>% 
    mutate(
      sim_idx = as.numeric(sub("Intervention_Patient_", "", name)) + 1
    )
  
  safe_max <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    return(max(x, na.rm = TRUE))
  }
  
  # 2. Flatten attributes (Joining by 'name' internally here is perfectly safe)
  sim_attributes <- get_mon_attributes(env) %>%
    group_by(name, key) %>% slice_tail(n = 1) %>% ungroup() %>%
    pivot_wider(names_from = key, values_from = value) %>%
    group_by(name) %>% 
    summarise(across(everything(), safe_max), .groups = "drop")
  
  # 3. Combine simulation outputs and keep our safe sim_idx
  sim_processed <- sim_arrivals %>%
    left_join(sim_attributes, by = "name") %>%
    select(sim_idx, outcome_status, call_time, appt_time, tier1_text_time, 
           tier2_text_time, wants_to_cancel, days_until_appt)
  
  # 4. Join back to the trial manifest using the deterministic sim_idx
  intervention_manifest <- trial_manifest %>% 
    filter(trial_arm == "Intervention") %>% 
    mutate(sim_idx = row_number()) %>% # Matches the +1 index we built above
    left_join(sim_processed, by = "sim_idx") %>% 
    select(-sim_idx) %>% # Clean up
    mutate(
      arrival_time = (week - 1) * mins_per_week,
      days_until_appt = ifelse(is.na(days_until_appt), runif(n(), min = 4, max = 14), days_until_appt),
      appt_time = ifelse(is.na(appt_time), arrival_time + (days_until_appt * mins_per_day), appt_time),
      tier1_text_time = ifelse(is.na(tier1_text_time), appt_time - (7 * mins_per_day), tier1_text_time),
      tier2_text_time = ifelse(is.na(tier2_text_time), appt_time - (3 * mins_per_day), tier2_text_time),
      outcome_status = replace_na(outcome_status, -1),
      sms_tier2_type = "Interactive Gate"
    )
  non_intervention_manifest <- trial_manifest %>%
    filter(trial_arm != "Intervention") %>%
    mutate(
      outcome_status = -2, call_time = NA_real_, days_until_appt = runif(n(), min = 4, max = 14),
      arrival_time = (week - 1) * mins_per_week, appt_time = arrival_time + (days_until_appt * mins_per_day),
      tier1_text_time = appt_time - (7 * mins_per_day), tier2_text_time = appt_time - (3 * mins_per_day),
      wants_to_cancel = 0, sms_tier2_type = "Passive Reminder"
    )
  
  # Step 4: Outcome Sampling & Compounding Risks
  final_trial_dataset <- bind_rows(intervention_manifest, non_intervention_manifest) %>%
    mutate(
      intervened_by_phone = ifelse(outcome_status == 1, TRUE, FALSE),
      is_truncated = ifelse(appt_time > total_sim_time, TRUE, FALSE),
      
      pathway_outcome = case_when(
        is_truncated ~ "Simulation Truncated",
        outcome_status == -1 ~ "Intervention Arm: Pending in Queue",
        trial_arm == "Not in Trial" ~ "Standard Low-Risk Track",
        trial_arm == "Control" ~ "Standard High-Risk Control Track",
        outcome_status == 0 ~ "Intervention Arm: Missed (Queue Timeout)",
        outcome_status == 1 & wants_to_cancel == 1 ~ "Intervention Arm: Cancelled via Coordinator",
        outcome_status == 1 & wants_to_cancel == 0 ~ "Intervention Arm: Confirmed via Coordinator"
      ),
      
      modulated_dna_prob = ml_baseline_risk * ifelse(appt_time > arrival_time, rr_tier1, 1) * ifelse(sms_tier2_type == "Passive Reminder", rr_tier2_passive, rr_tier2_interact) * ifelse(pathway_outcome == "Intervention Arm: Confirmed via Coordinator", rr_tier3_call, 1),
      
      modulated_dna_prob = pmin(1, pmax(0, modulated_dna_prob)),
      
      final_attendance = case_when(
        pathway_outcome == "Simulation Truncated" ~ "Censored/Right-Truncated",
        pathway_outcome == "Intervention Arm: Pending in Queue" ~ "Censored/Right-Truncated",
        pathway_outcome == "Intervention Arm: Cancelled via Coordinator" ~ "Cancelled/Rearranged",
        runif(n()) < modulated_dna_prob ~ "DNA",
        TRUE ~ "Attended"
      )
    )
  
  return(final_trial_dataset)
}

plot_trial_journeys <- function(trial_data, num_samples_per_group = 3) {
  
  # 1. Sample patients and set up ordering
  sampled_patients <- trial_data %>%
    group_by(trial_arm, final_attendance) %>%
    slice_sample(n = num_samples_per_group, replace = TRUE) %>%
    ungroup() %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    mutate(patient_id = factor(patient_id)) %>%
    mutate(patient_plot_order = reorder(patient_id, appt_time))
  
  # 2. Extract milestone points for scatter layer
  touchpoint <- sampled_patients %>%
    transmute(
      patient_id = patient_plot_order,
      trial_arm = trial_arm,
      final_state = final_attendance,
      List_Extraction = arrival_time / 60,
      Tier1_Automated_Text = tier1_text_time / 60,
      Tier2_Passive_Reminder = ifelse(sms_tier2_type == "Passive Reminder", tier2_text_time / 60, NA_real_),
      Tier2_Interactive_Gate = ifelse(sms_tier2_type == "Interactive Gate", tier2_text_time / 60, NA_real_),
      Coordinator_Call = call_time / 60,
      Appointment = appt_time / 60
    ) %>%
    pivot_longer(
      cols = c(List_Extraction, Tier1_Automated_Text, Tier2_Passive_Reminder, Tier2_Interactive_Gate, Coordinator_Call, Appointment),
      names_to = "Milestone",
      values_to = "Time_Hours"
    ) %>%
    filter(!is.na(Time_Hours)) %>%
    filter(Time_Hours >= 0) # Only display milestone points within the study window
  
  # 3. Extract journey lifespans for segment layer
  lifespans <- sampled_patients %>%
    transmute(
      patient_id = patient_plot_order,
      trial_arm = trial_arm,
      final_state = final_attendance,
      start = pmax(0, pmin(arrival_time, tier1_text_time, na.rm = TRUE) / 60),
      end = appt_time / 60
    )
  
  # 4. Build plot
  p <- ggplot() +
    # Journey duration lines
    geom_segment(
      data = lifespans,
      aes(x = start, xend = end, y = patient_id, yend = patient_id, color = final_state),
      linewidth = 1.2, alpha = 0.75
    ) +
    # Event milestone markers
    geom_point(
      data = touchpoint,
      alpha = 0.5,
      aes(x = Time_Hours, y = patient_id, shape = Milestone, fill = Milestone),
      size = 3.2, stroke = 0.6, color = "black"
    ) +
    # Custom shapes (shapes 21-25 allow separate fill and border colors)
    scale_shape_manual(values = c(
      "List_Extraction"        = 21, # Circle
      "Tier1_Automated_Text"   = 22, # Square
      "Tier2_Passive_Reminder" = 23, # Diamond
      "Tier2_Interactive_Gate" = 24, # Triangle Up
      "Coordinator_Call"       = 25, # Triangle Down
      "Appointment"            = 21  # Circle
    )) +
    # Distinct colors for touchpoint
    scale_fill_manual(values = c(
      "List_Extraction"        = "#66c2a5",
      "Tier1_Automated_Text"   = "#fc8d62",
      "Tier2_Passive_Reminder" = "#8da0cb",
      "Tier2_Interactive_Gate" = "#e78ac3",
      "Coordinator_Call"       = "#a6d854",
      "Appointment"            = "#ffd92f"
    )) +
    scale_colour_manual(values = c(
      "Attended" = "darkolivegreen3", 
      "Cancelled/Rearranged" = "cornsilk4",
      "DNA" = "coral2",
      "Censored/Right-Truncated" = "darkslategray4"
      )) +
    # Facet by trial arm so comparisons are easy
    facet_wrap(~ trial_arm, scales = "free_y", ncol = 1) +
    # Layout and labels
    labs(
      title = "Patient Journey Timelines across Trial Arms",
      x = "Time Relative to Appointment (Hours)",
      y = "Patient ID",
      color = "Final Attendance",
      shape = "Touchpoint",
      fill = "Touchpoint"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
  
  return(p)
}

# Execute clean data population
trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = 12)
plot_trial_journeys(trial_data)



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

calculate_simulation_power <- function(target_weeks, iterations = 100, target_effect = 0.40) {
  significant_results <- 0
  
  for (i in 1:iterations) {
    # 1. Run the simulation
    sim_data <- simulate_clinical_trial_advanced(
      weeks_to_simulate = target_weeks,
      coordinator_capacity = 1 
    )
    
    # 2. Filter and apply the New Factor Structure
    analysis_data <- sim_data %>%
      filter(!final_attendance %in% c("Censored/Right-Truncated", "Cancelled/Rearranged")) %>%
      filter(trial_arm != "Not in Trial") %>%
      mutate(
        dna_outcome = ifelse(final_attendance == "DNA", 1, 0),
        pp_exposure = case_when(
          trial_arm == "Control" ~ "1_Control", 
          trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "2_Tier2_Interactive_Only",
          trial_arm == "Intervention" & intervened_by_phone == TRUE ~ "3_Tier3_Call_Reached"
        ),
        pp_exposure = factor(pp_exposure, levels = c("1_Control", "2_Tier2_Interactive_Only", "3_Tier3_Call_Reached"))
      )
    
    # 3. Fit Modified Poisson Model
    fit <- glm(dna_outcome ~ pp_exposure + ml_baseline_risk,
               data = analysis_data, family = poisson(link = "log"))
    
    # 4. Extract Robust P-Values
    robust_se <- coeftest(fit, vcov = sandwich)
    tidy_results <- tidy(robust_se)
    
    # 5. Check if target effect was detected successfully (UPDATED NAME)
    call_p_value <- tidy_results %>% 
      filter(term == "pp_exposure3_Tier3_Call_Reached") %>% 
      pull(p.value)
    
    if (length(call_p_value) > 0 && !is.na(call_p_value) && call_p_value < 0.05) {
      significant_results <- significant_results + 1
    }
  }
  
  power <- significant_results / iterations
  return(power)
}

# Example Usage: (Run this with varying target_weeks to find when power >= 0.80)
estimated_power <- calculate_simulation_power(target_weeks = 12, iterations = 50)
print(paste("Power at 12 weeks:", estimated_power))


plot_recovered_parameters <- function(pp_model_results) {
  
  # Extract estimates and robust standard errors from the coeftest object
  plot_data <- tidy(pp_model_results) %>%
    # UPDATED: Filter for the new factor coefficient names
    filter(term %in% c("pp_exposure2_Tier2_Interactive_Only", "pp_exposure3_Tier3_Call_Reached")) %>%
    mutate(
      Risk_Ratio = exp(estimate),
      Conf_Low = exp(estimate - 1.96 * std.error),
      Conf_High = exp(estimate + 1.96 * std.error),
      
      # UPDATED: Map to the new names
      True_Parameter = case_when(
        term == "pp_exposure2_Tier2_Interactive_Only" ~ 0.85,
        term == "pp_exposure3_Tier3_Call_Reached" ~ 0.40
      ),
      Term_Label = case_when(
        term == "pp_exposure2_Tier2_Interactive_Only" ~ "Tier 2: Interactive Text Gate\n(True RR = 0.85)",
        term == "pp_exposure3_Tier3_Call_Reached" ~ "Tier 3: Coordinator Call\n(True RR = 0.40)"
      )
    )
  
  # Generate Forest Plot (Unchanged)
  ggplot(plot_data, aes(x = Risk_Ratio, y = Term_Label)) +
    geom_pointrange(aes(xmin = Conf_Low, xmax = Conf_High), 
                    color = "#2c3e50", size = 1, linewidth = 1.2) +
    geom_point(aes(x = True_Parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    scale_x_continuous(limits = c(0.2, 1.2), breaks = seq(0.2, 1.2, 0.2)) +
    labs(
      title = "Validation of Pathway Component Effect Sizes",
      subtitle = "Red diamonds represent programmed ground truth. Black ranges are recovered 95% CIs.",
      x = "Recovered Risk Ratio (RR)",
      y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

# Execute Plot
plot_recovered_parameters(pp_results)
