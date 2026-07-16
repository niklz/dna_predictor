library(simmer)
library(dplyr)
library(tidyr)

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
    rr_tier1          = 0.95,  # 5% relative reduction from the universal text
    rr_tier2_passive  = 0.90,  # 10% relative reduction from passive text
    rr_tier2_interact = 0.85,  # 15% relative reduction from interactive text setup
    rr_tier3_call     = 0.40   # 60% relative reduction if they talk to a human
) {
  
  mins_per_day <- 1440
  mins_per_week <- 7 * mins_per_day
  total_sim_time <- weeks_to_simulate * mins_per_week
  work_mins_per_week <- 5 * 8 * 60
  mean_handling_time <- work_mins_per_week / max_calls_per_week
  
  # ---------------------------------------------------------
  # STEP 1: Global Trial Cohort & ML Risk Generation
  # ---------------------------------------------------------
  num_weeks <- weeks_to_simulate
  weekly_total_sizes <- rpois(num_weeks, lambda = total_patients_per_week)
  
  trial_manifest <- data.frame(week = rep(1:num_weeks, times = weekly_total_sizes)) %>%
    mutate(
      patient_id = paste0("PT_", row_number()),
      risk_profile = ifelse(runif(n()) < ml_high_risk_prop, "High-Risk", "Low-Risk"),
      
      # NEW: Direct ML Model Risk Output (Continuous known baseline probability)
      # Low-risk hover around 4%; High-risk hover around 20% baseline
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
  
  # ---------------------------------------------------------
  # STEP 2: Run SIMMER for Intervention Arm
  # ---------------------------------------------------------
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
  
  # ---------------------------------------------------------
  # STEP 3: Extract & Flatten Sim Data (SAFE FROM WARNINGS)
  # ---------------------------------------------------------
  sim_arrivals <- get_mon_arrivals(env) %>% mutate(join_idx = row_number())
  
  # A custom helper to handle completely missing columns without throwing -Inf warnings
  safe_max <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    return(max(x, na.rm = TRUE))
  }
  
  sim_attributes <- get_mon_attributes(env) %>%
    group_by(name, key) %>% slice_tail(n = 1) %>% ungroup() %>%
    pivot_wider(names_from = key, values_from = value) %>%
    group_by(name) %>% 
    # Use our safe helper function instead of raw max()
    summarise(across(everything(), safe_max), .groups = "drop")
  
  sim_processed <- sim_arrivals %>%
    left_join(sim_attributes, by = "name") %>%
    arrange(arrival_time) %>% mutate(join_idx = row_number()) %>%
    select(join_idx, outcome_status, call_time, appt_time, tier1_text_time, tier2_text_time, wants_to_cancel, days_until_appt)
  
  intervention_manifest <- trial_manifest %>% 
    filter(trial_arm == "Intervention") %>% mutate(join_idx = row_number()) %>%
    left_join(sim_processed, by = "join_idx") %>% select(-join_idx) %>%
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
  
  # ---------------------------------------------------------
  # STEP 4: Robust Incremental Risk Calculation & Sampling
  # ---------------------------------------------------------
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
      
      # COMPUTE PATIENT-SPECIFIC MODULATED RISK
      # Start with baseline risk, multiply by each achieved milestone's Risk Ratio
      modulated_dna_prob = ml_baseline_risk * # Tier 1 applied universally if appointment wasn't immediately truncated
        ifelse(appt_time > arrival_time, rr_tier1, 1) * # Tier 2 variation
        ifelse(sms_tier2_type == "Passive Reminder", rr_tier2_passive, rr_tier2_interact) * # Tier 3 human intervention multiplier
        ifelse(pathway_outcome == "Intervention Arm: Confirmed via Coordinator", rr_tier3_call, 1),
      
      # Boundary protection to ensure probability stays between 0 and 1
      modulated_dna_prob = pmin(1, pmax(0, modulated_dna_prob)),
      
      # RANDOM OUTCOME SAMPLING BASED ON INDIVIDUAL RISK
      final_attendance = case_when(
        pathway_outcome == "Simulation Truncated" ~ "Censored/Right-Truncated",
        pathway_outcome == "Intervention Arm: Pending in Queue" ~ "Censored/Right-Truncated",
        pathway_outcome == "Intervention Arm: Cancelled via Coordinator" ~ "Cancelled/Rearranged",
        
        # Draw a Bernoulli trial using each patient's custom calculated probability
        runif(n()) < modulated_dna_prob ~ "DNA",
        TRUE ~ "Attended"
      )
    )
  
  return(final_trial_dataset)
}


plot_trial_journeys <- function(trial_data, num_samples_per_group = 3) {
  
  sampled_patients <- trial_data %>%
    group_by(trial_arm, final_attendance) %>%
    slice_sample(n = num_samples_per_group, replace = TRUE) %>% 
    ungroup() %>%
    distinct(patient_id, .keep_all = TRUE) %>% 
    mutate(patient_id = factor(patient_id)) %>%
    mutate(patient_plot_order = reorder(patient_id, appt_time))
  
  milestones <- sampled_patients %>%
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
  
  lifespans <- sampled_patients %>%
    transmute(
      patient_id = patient_plot_order,
      trial_arm = trial_arm,
      final_state = final_attendance,
      # CRITICAL BUG FIX: Use pmax(0, ...) so the line clips to the edge instead of breaking
      start = pmax(0, pmin(arrival_time, tier1_text_time) / 60), 
      end = appt_time / 60
    )
  
  ggplot() +
    geom_segment(
      data = lifespans,
      aes(x = start, xend = end, y = patient_id, yend = patient_id, color = final_state),
      linewidth = 1.2, alpha = 0.5
    ) +
    geom_point(
      data = milestones,
      aes(x = Time_Hours, y = patient_id, shape = Milestone, fill = Milestone),
      size = 3.2, stroke = 0.6, color = "black"
    ) +
    scale_shape_manual(values = c(
      "List_Extraction" = 21, 
      "Tier1_Automated_Text" = 22, 
      "Tier2_Passive_Reminder" = 23,  
      "Tier2_Interactive_Gate" = 23,  
      "Coordinator_Call" = 25, 
      "Appointment" = 24
    )) +
    scale_fill_manual(values = c(
      "List_Extraction" = "#999999", 
      "Tier1_Automated_Text" = "#CC79A7", 
      "Tier2_Passive_Reminder" = "#56B4E9", 
      "Tier2_Interactive_Gate" = "#0072B2", 
      "Coordinator_Call" = "#E69F00", 
      "Appointment" = "#D55E00"
    )) +
    scale_color_manual(values = c(
      "Attended" = "#009E73",         
      "DNA" = "#D55E00",              
      "Cancelled/Rearranged" = "#0072B2",
      "Censored/Right-Truncated" = "#999999"
    )) +
    facet_wrap(~ trial_arm, scales = "free_y", nrow = 1) +
    labs(
      title = "Synthetic Trial Dataset: Evaluated Patient Cohorts",
      subtitle = "Fixed segment clipping for historical text message boundaries.",
      x = "Simulation Timeline (Hours)",
      y = "Sampled Patients per Cohort",
      color = "Ground-Truth Final Status",
      shape = "Timeline Milestones",
      fill = "Timeline Milestones"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      strip.text = element_text(face = "bold", size = 11)
    )
}


# Change this line at the bottom of your script:
trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = 12)

# Then run the plot function right after
plot_trial_journeys(trial_data, num_samples_per_group = 8)


