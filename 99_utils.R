library(ggplot2)
library(ggalluvial)
library(broom)
library(dplyr)
library(purrr)

simulate_clinical_trial_advanced <- function(
    weeks_to_simulate = 8,
    total_patients_per_week = 1000,
    ml_high_risk_prop = 0.25,
    trial_allocation_ratio = 0.50,
    
    # Coordinator Constraints (Applies ONLY to High-Risk Intervention)
    coordinator_capacity = 1,
    max_calls_per_week = 100,
    response_prob = 0.30,
    prob_wish_to_cancel = 0.20,
    
    # INCREMENTAL EFFECT SIZES (Risk Ratios)
    rr_tier1          = 0.95,  
    rr_tier2_passive  = 0.90,  
    rr_tier2_interact = 0.85,  
    rr_tier3_call     = 0.40   
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
  
  # Split cleanly into Intervention vs Non-Intervention (Control + Low-Risk)
  intervention_df <- trial_manifest %>% filter(trial_arm == "Intervention") %>% mutate(sim_idx = row_number())
  standard_df <- trial_manifest %>% filter(trial_arm != "Intervention") %>% mutate(sim_idx = row_number())
  
  intervention_arrivals <- (intervention_df$week - 1) * mins_per_week
  standard_arrivals <- (standard_df$week - 1) * mins_per_week
  
  # Step 2: DES Simulation Process
  env <- simmer("trial_pathway")
  
  intervention_path <- trajectory("High_Risk_Intervention_Journey") %>%
    set_attribute("arrival_time", function() { simmer::now(env) }) %>%
    set_attribute("days_until_appt", function() { runif(1, min = 4, max = 14) }) %>%
    set_attribute("appt_time", function() { simmer::now(env) + (get_attribute(env, "days_until_appt") * mins_per_day) }) %>%
    
    # FIX 1: Prevent Tier 1 from time-traveling before trial entry
    set_attribute("tier1_text_time", function() { 
      base_t1 <- get_attribute(env, "appt_time") - (7 * mins_per_day) + rnorm(1, 0, 120)
      max(simmer::now(env), base_t1) 
    }) %>%
    
    # FIX 2: Physically advance the clock before sending Tier 2
    timeout(function() { 
      target_t2 <- get_attribute(env, "appt_time") - (3 * mins_per_day) + rnorm(1, 0, 60)
      max(0, target_t2 - simmer::now(env)) 
    }) %>%
    set_attribute("tier2_text_time", function() { simmer::now(env) }) %>%
    set_attribute("response_delay_hours", function() { rexp(1, rate = 1/12) }) %>% 
    
    branch(
      function() {
        appt_t <- get_attribute(env, "appt_time")
        curr_t <- simmer::now(env)
        delay_mins <- get_attribute(env, "response_delay_hours") * 60
        response_t <- curr_t + delay_mins
        cutoff_t <- appt_t - mins_per_day 
        
        did_respond <- rbinom(1, 1, response_prob)
        
        if (did_respond == 0 || response_t > cutoff_t) { return(3) }
        return(sample(c(1, 2), size = 1, prob = c(0.70, 0.30))) 
      },
      continue = c(TRUE, TRUE, TRUE),
      
      trajectory("Intent_Confirm") %>%
        set_attribute("response_type", 1) %>%
        set_attribute("outcome_status", 2),
      
      trajectory("Intent_Reschedule_Cancel") %>%
        set_attribute("response_type", 2) %>%
        # FIX 3a: Wait out the response delay before hitting the queue
        timeout(function() { get_attribute(env, "response_delay_hours") * 60 }) %>%
        renege_in(
          t = function() { max(0, get_attribute(env, "appt_time") - simmer::now(env)) },
          out = trajectory() %>% set_attribute("outcome_status", 0) 
        ) %>%
        seize("coordinator", 1) %>%
        renege_abort() %>% 
        set_attribute("call_time", function() { simmer::now(env) }) %>%
        timeout(function() { max(5, rnorm(1, mean = mean_handling_time, sd = 5)) }) %>%
        set_attribute("wants_to_cancel", 1) %>% 
        release("coordinator", 1) %>%
        set_attribute("outcome_status", 1),
      
      trajectory("Intent_No_Response") %>%
        set_attribute("response_type", 3) %>%
        # FIX 3b: Wait until exactly the 24-hour cutoff before triggering chase
        timeout(function() { 
          max(0, get_attribute(env, "appt_time") - mins_per_day - simmer::now(env)) 
        }) %>%
        renege_in(
          t = function() { max(0, get_attribute(env, "appt_time") - simmer::now(env)) },
          out = trajectory() %>% set_attribute("outcome_status", 0) 
        ) %>%
        seize("coordinator", 1) %>%
        renege_abort() %>% 
        set_attribute("call_time", function() { simmer::now(env) }) %>%
        timeout(function() { max(5, rnorm(1, mean = mean_handling_time, sd = 5)) }) %>%
        set_attribute("wants_to_cancel", function() { rbinom(1, 1, prob_wish_to_cancel) }) %>%
        release("coordinator", 1) %>%
        set_attribute("outcome_status", 1)
    )
  
  # Ensure standard paths adhere to chronological limits as well
  standard_path <- trajectory("Standard_Care_Journey") %>%
    set_attribute("arrival_time", function() { simmer::now(env) }) %>%
    set_attribute("days_until_appt", function() { runif(1, min = 4, max = 14) }) %>%
    set_attribute("appt_time", function() { simmer::now(env) + (get_attribute(env, "days_until_appt") * mins_per_day) }) %>%
    set_attribute("tier1_text_time", function() { 
      base_t1 <- get_attribute(env, "appt_time") - (7 * mins_per_day) + rnorm(1, 0, 120)
      max(simmer::now(env), base_t1)
    }) %>%
    set_attribute("tier2_text_time", function() { 
      base_t2 <- get_attribute(env, "appt_time") - (3 * mins_per_day) + rnorm(1, 0, 60)
      max(get_attribute(env, "tier1_text_time"), base_t2) # Must happen after Tier 1
    })
  
  working_hours <- schedule(timetable = c(0, 9*60, 17*60), values = c(0, coordinator_capacity, 0), period = mins_per_day)
  
  env %>%
    add_resource("coordinator", capacity = working_hours) %>%
    add_generator("Intervention_", intervention_path, at(intervention_arrivals), mon = 2) %>%
    add_generator("Standard_", standard_path, at(standard_arrivals), mon = 2) %>%
    run(until = total_sim_time)
  
  sim_arrivals <- get_mon_arrivals(env) %>% 
    mutate(
      arm_group = sub("_.*", "", name),
      sim_idx = as.numeric(sub(".*_", "", name)) + 1 
    )
  
  safe_max <- function(x) { if (all(is.na(x))) return(NA_real_) else return(max(x, na.rm = TRUE)) }
  
  sim_attributes <- get_mon_attributes(env) %>%
    group_by(name, key) %>% slice_tail(n = 1) %>% ungroup() %>%
    pivot_wider(names_from = key, values_from = value) %>%
    group_by(name) %>% summarise(across(everything(), safe_max), .groups = "drop")
  
  sim_processed <- sim_arrivals %>% left_join(sim_attributes, by = "name")
  
  expected_cols <- c("outcome_status", "call_time", "appt_time", "tier1_text_time", "tier2_text_time", "wants_to_cancel", "days_until_appt")
  for (col in expected_cols) { if (!col %in% names(sim_processed)) sim_processed[[col]] <- NA_real_ }
  
  intervention_manifest <- intervention_df %>%
    left_join(sim_processed %>% filter(arm_group == "Intervention"), by = "sim_idx") %>% 
    dplyr::select(-sim_idx, -name, -arm_group) %>% 
    mutate(
      outcome_status = replace_na(outcome_status, -1),
      sms_tier2_type = "Interactive Gate"
    )
  
  non_intervention_manifest <- standard_df %>%
    left_join(sim_processed %>% filter(arm_group == "Standard"), by = "sim_idx") %>%
    dplyr::select(-sim_idx, -name, -arm_group) %>% 
    mutate(
      outcome_status = -2, 
      wants_to_cancel = 0, 
      sms_tier2_type = "Passive Reminder",
      call_time = NA_real_ 
    )
  
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
        outcome_status == 1 & wants_to_cancel == 0 ~ "Intervention Arm: Confirmed via Coordinator",
        TRUE ~ "Completed Auto-Text Only"
      ),
      
      modulated_dna_prob = ml_baseline_risk * ifelse(appt_time > arrival_time, rr_tier1, 1) * ifelse(sms_tier2_type == "Passive Reminder", rr_tier2_passive, rr_tier2_interact) * ifelse(pathway_outcome == "Intervention Arm: Confirmed via Coordinator", rr_tier3_call, 1),
      
      modulated_dna_prob = pmin(1, pmax(0, modulated_dna_prob)),
      
      final_attendance = case_when(
        pathway_outcome %in% c("Simulation Truncated", "Intervention Arm: Pending in Queue", "Intervention Arm: Missed (Queue Timeout)") ~ "Censored",
        pathway_outcome == "Intervention Arm: Cancelled via Coordinator" ~ "Cancelled",
        runif(n()) < modulated_dna_prob ~ "DNA",
        TRUE ~ "Attended"
      )
    )
  
  # Final Dataset - Preserving all columns required for mapping
  sop_dataset <- final_trial_dataset %>%
    dplyr::select(
      patient_id, trial_arm, risk_profile, ml_baseline_risk,       
      sms_tier2_type, intervened_by_phone, final_attendance,       
      arrival_time, appt_time, tier1_text_time, tier2_text_time, call_time, pathway_outcome 
    )
  
  return(sop_dataset)
}

plot_trial_journeys <- function(trial_data, num_samples_per_group = 3) {
  
  sampled_patients <- trial_data %>%
    filter(!final_attendance %in% c("Censored", "Cancelled")) %>%
    group_by(trial_arm, final_attendance) %>%
    slice_sample(n = num_samples_per_group, replace = TRUE) %>%
    ungroup() %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    mutate(patient_plot_order = reorder(as.factor(patient_id), appt_time))
  
  touchpoint <- sampled_patients %>%
    transmute(
      patient_id = patient_plot_order, trial_arm,
      List_Extraction = arrival_time / 60,
      Tier1_Automated_Text = tier1_text_time / 60,
      Tier2_Passive_Reminder = ifelse(sms_tier2_type == "Passive Reminder", tier2_text_time / 60, NA_real_),
      Tier2_Interactive_Gate = ifelse(sms_tier2_type == "Interactive Gate", tier2_text_time / 60, NA_real_),
      Coordinator_Call = call_time / 60, Appointment = appt_time / 60
    ) %>%
    pivot_longer(cols = -c(patient_id, trial_arm), names_to = "Milestone", values_to = "Time_Hours") %>%
    filter(!is.na(Time_Hours) & Time_Hours >= 0)
  
  lifespans <- sampled_patients %>%
    transmute(patient_id = patient_plot_order, trial_arm, final_attendance,
              start = pmax(0, pmin(arrival_time, tier1_text_time, na.rm = TRUE) / 60),
              end = appt_time / 60)
  
  max_hrs <- max(lifespans$end, na.rm = TRUE)
  
  weekend_starts <- seq(120, max_hrs, by = 168)
  weekends <- data.frame(
    xmin = weekend_starts,
    xmax = weekend_starts + 48
  ) %>% filter(xmin <= max_hrs)
  
  ggplot() +
    geom_rect(data = weekends, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey90", alpha = 0.6, inherit.aes = FALSE) +
    geom_segment(data = lifespans, aes(x = start, xend = end, y = patient_id, yend = patient_id, color = final_attendance),
                 linewidth = 1.2, alpha = 0.75) +
    geom_point(data = touchpoint, aes(x = Time_Hours, y = patient_id, shape = Milestone, fill = Milestone),
               size = 3.2, stroke = 0.6, color = "black", alpha = 0.8) +
    scale_shape_manual(values = c("List_Extraction"=21, "Tier1_Automated_Text"=22, "Tier2_Passive_Reminder"=23, 
                                  "Tier2_Interactive_Gate"=24, "Coordinator_Call"=25, "Appointment"=21)) +
    scale_fill_manual(values = c("List_Extraction"="#66c2a5", "Tier1_Automated_Text"="#fc8d62", 
                                 "Tier2_Passive_Reminder"="#8da0cb", "Tier2_Interactive_Gate"="#e78ac3", 
                                 "Coordinator_Call"="#a6d854", "Appointment"="#ffd92f")) +
    scale_colour_manual(values = c("Attended" = "darkolivegreen3", "DNA" = "coral2")) +
    facet_wrap(~ trial_arm, scales = "free_y", ncol = 1) +
    labs(title = "Patient Journey Timelines", subtitle = "Grey bands indicate weekends.",
         x = "Time (Hours from simulation start)", y = "") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", axis.text.y = element_blank())
}


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


plot_recovered_parameters <- function(
    pp_model_results, 
    rr_tier2_passive  = 0.90, 
    rr_tier2_interact = 0.85, 
    rr_tier3_call     = 0.40,
    rr_tier1          = 0.95
) {
  
  # Calculate the true relative ratios expected by the GLM offset design:
  # - Tier 2 is relative to the Control group (Passive Reminder)
  # - Tier 3 cumulative effect includes Tier 1 * Tier 2 Interact * Tier 3 Call, relative to Control (Tier 1 * Tier 2 Passive)
  true_rr_tier2 <- rr_tier2_interact / rr_tier2_passive
  true_rr_tier3 <- (rr_tier1 * rr_tier2_interact * rr_tier3_call) / (rr_tier1 * rr_tier2_passive)
  
  plot_data <- tidy(pp_model_results) %>%
    filter(term %in% c("pp_exposure2_Tier2_Interactive_Only", "pp_exposure3_Tier3_Call_Reached")) %>%
    mutate(
      Risk_Ratio = exp(estimate),
      Conf_Low = exp(estimate - 1.96 * std.error),
      Conf_High = exp(estimate + 1.96 * std.error),
      
      Tier = case_when(
        grepl("Tier2", term) ~ "Tier 2 Impacts",
        grepl("Tier3", term) ~ "Tier 3 Impacts"
      ),
      True_Parameter = case_when(
        grepl("Tier2", term) ~ true_rr_tier2,
        grepl("Tier3", term) ~ true_rr_tier3
      ),
      Term_Label = case_when(
        grepl("Tier2", term) ~ paste0("Interactive gate vs passive reminder\n(true RR = ", round(true_rr_tier2, 3), ")"),
        grepl("Tier3", term) ~ paste0("Coordinator call vs text only\n(true RR = ", round(true_rr_tier3, 3), ")")
      )
    )
  
  ggplot(plot_data, aes(x = Risk_Ratio, y = Term_Label)) +
    geom_pointrange(aes(xmin = Conf_Low, xmax = Conf_High), 
                    color = "#2c3e50", size = 1, linewidth = 1.2) +
    geom_point(aes(x = True_Parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    facet_grid(Tier ~ ., scales = "free_y", space = "free_y") +
    # scale_x_continuous(limits = c(0.2, 1.2), breaks = seq(0.2, 1.2, 0.2)) +
    labs(
      title = "Validation of pathway component effect sizes",
      subtitle = "Red diamonds represent programmed relative ground truth. Black ranges are recovered 95% CIs.",
      x = "Recovered risk ratio (RR)", y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = NA),
      plot.title = element_text(face = "bold")
    )
}

# =========================================================================
# SECTION 4: PROCESS FLOW & BRANCH VOLUME METRICS
# =========================================================================
analyze_pathway_volumes <- function(trial_data) {
  volume_summary <- trial_data %>%
    count(trial_arm, risk_profile, sms_tier2_type, intervened_by_phone, final_attendance) %>%
    mutate(percentage = round(n / sum(n) * 100, 2))
  
  print("--- PATIENT FLOW VOLUME BREAKDOWN ---")
  print(volume_summary)
  return(volume_summary)
}

plot_transition_time_distributions <- function(trial_data) {
  
  transition_data <- trial_data %>%
    filter(!final_attendance %in% c("Censored", "Cancelled")) %>%
    mutate(
      # Calculate true forward operational latencies (in hours)
      dur_entry_to_t1 = (tier1_text_time - arrival_time) / 60,
      dur_t1_to_t2    = (tier2_text_time - tier1_text_time) / 60,
      dur_t2_to_call  = ifelse(intervened_by_phone == TRUE, (call_time - tier2_text_time) / 60, NA_real_),
      dur_call_to_appt = ifelse(intervened_by_phone == TRUE, (appt_time - call_time) / 60, (appt_time - tier2_text_time) / 60)
    ) %>%
    dplyr::select(patient_id, trial_arm, dur_entry_to_t1, dur_t1_to_t2, dur_t2_to_call, dur_call_to_appt) %>%
    pivot_longer(
      cols = starts_with("dur_"),
      names_to = "Transition_Stage",
      values_to = "Hours"
    ) %>%
    filter(!is.na(Hours) & Hours >= 0) %>%
    mutate(
      Stage_Label = case_when(
        Transition_Stage == "dur_entry_to_t1"  ~ "1. Trial Entry -> Tier 1 Text",
        Transition_Stage == "dur_t1_to_t2"     ~ "2. Tier 1 Text -> Tier 2 Touchpoint",
        Transition_Stage == "dur_t2_to_call"   ~ "3. Tier 2 Touchpoint -> Coordinator Call",
        Transition_Stage == "dur_call_to_appt" ~ "4. Final Action -> Appointment"
      )
    )
  
  ggplot(transition_data, aes(x = Hours, fill = trial_arm)) +
    geom_density(alpha = 0.5, color = NA, adjust = 1.5) +
    facet_wrap(~ Stage_Label, scales = "free") +
    labs(
      title = "Distribution of Transition Times Between Pathway Stages",
      subtitle = "Evaluates operational latency and variance across trial arms",
      x = "Duration (Hours)",
      y = "Density",
      fill = "Trial Arm"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.text = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 11)
    )
}


plot_pathway_gantt_chart <- function(trial_data, max_patients_to_plot = 150) {
  
  # Ensure optional columns exist safely to prevent missing object errors
  if (!"wants_to_cancel" %in% names(trial_data)) {
    trial_data$wants_to_cancel <- 0
  }
  
  # Prepare and filter dataset for plotting
  plot_data_prep <- trial_data %>%
    filter(final_attendance != "Censored") %>%
    mutate(
      branch_category = case_when(
        trial_arm == "Control" ~ "Control Arm (Standard Track)",
        trial_arm == "Not in Trial" ~ "Low-Risk Standard Track",
        trial_arm == "Intervention" & intervened_by_phone == TRUE & wants_to_cancel == 1 ~ "Intervention: Cancelled via Call",
        trial_arm == "Intervention" & intervened_by_phone == TRUE & wants_to_cancel == 0 ~ "Intervention: Confirmed via Call",
        trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "Intervention: Coordinator Queue Timeout / Text-Only",
        TRUE ~ "Other Path"
      )
    ) %>%
    group_by(branch_category) %>%
    slice_head(n = max_patients_to_plot) %>% 
    ungroup() %>%
    mutate(patient_plot_order = reorder(as.factor(patient_id), appt_time))
  
  # Extract milestone touchpoints in hours from simulation start
  touchpoints <- plot_data_prep %>%
    transmute(
      patient_id = patient_plot_order,
      branch_category,
      final_attendance,
      List_Extraction = arrival_time / 60,
      Tier1_Text = tier1_text_time / 60,
      Tier2_Touchpoint = tier2_text_time / 60,
      Coordinator_Call = call_time / 60,
      Appointment = appt_time / 60
    ) %>%
    pivot_longer(cols = -c(patient_id, branch_category, final_attendance), 
                 names_to = "Milestone", values_to = "Time_Hours") %>%
    filter(!is.na(Time_Hours) & Time_Hours >= 0)
  
  # Calculate patient lifespan bars
  lifespans <- plot_data_prep %>%
    transmute(
      patient_id = patient_plot_order,
      branch_category,
      final_attendance,
      start = pmax(0, pmin(arrival_time, tier1_text_time, na.rm = TRUE) / 60),
      end = appt_time / 60
    )
  
  # Weekend background blocks calculation (Safe length alignment)
  max_hrs <- max(lifespans$end, na.rm = TRUE)
  weekend_starts <- seq(120, max_hrs, by = 168)
  
  weekends <- data.frame(
    xmin = weekend_starts,
    xmax = pmin(weekend_starts + 48, max_hrs + 168) # Ensures xmax never mismatches xmin length
  ) %>% filter(xmin <= max_hrs)
  
  # Build the Gantt-style visualization
  ggplot() +
    geom_rect(data = weekends, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey95", alpha = 0.8, inherit.aes = FALSE) +
    geom_segment(data = lifespans, 
                 aes(x = start, xend = end, y = patient_id, yend = patient_id, color = final_attendance),
                 linewidth = 1.0, alpha = 0.6) +
    geom_point(data = touchpoints, 
               aes(x = Time_Hours, y = patient_id, shape = Milestone, fill = Milestone),
               size = 2.2, stroke = 0.4, color = "black", alpha = 0.85) +
    scale_shape_manual(values = c(
      "List_Extraction" = 21, 
      "Tier1_Text" = 22, 
      "Tier2_Touchpoint" = 23, 
      "Coordinator_Call" = 25, 
      "Appointment" = 21
    )) +
    scale_fill_manual(values = c(
      "List_Extraction" = "#66c2a5", 
      "Tier1_Text" = "#fc8d62", 
      "Tier2_Touchpoint" = "#8da0cb", 
      "Coordinator_Call" = "#e78ac3", 
      "Appointment" = "#ffd92f"
    )) +
    scale_colour_manual(values = c("Attended" = "#2ecc71", "DNA" = "#e74c3c", "Cancelled" = "#95a5a6")) +
    facet_grid(branch_category ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "Patient Pathway Gantt Chart: Branch Volume Breakdown",
      subtitle = "Timelines trace individual patient journeys down distinct interactive text and coordinator queues. Grey bars indicate weekends.",
      x = "Simulation Time (Hours)",
      y = "Patients (Grouped by Final Branch)",
      color = "Clinical Outcome",
      fill = "Milestone",
      shape = "Milestone"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey90", color = NA),
      axis.text.y = element_blank(),
      legend.position = "bottom"
    )
}

plot_pathway_flow <- function(trial_data) {
  
  # Ensure optional columns exist safely
  if (!"wants_to_cancel" %in% names(trial_data)) {
    trial_data$wants_to_cancel <- 0
  }
  
  # 1. Aggregate and categorize the pathway branches
  flow_data <- trial_data %>%
    filter(final_attendance != "Censored") %>%
    mutate(
      # Stage 1: Initial Arm Assignment
      Stage1_Arm = case_when(
        trial_arm == "Control" ~ "Control (High Risk)",
        trial_arm == "Not in Trial" ~ "Standard (Low Risk)",
        TRUE ~ "Intervention"
      ),
      
      # Stage 2: Interaction / Queue Outcome
      Stage2_Action = case_when(
        Stage1_Arm != "Intervention" ~ "Standard Track",
        intervened_by_phone == FALSE ~ "Resolved via Text / Timeout",
        intervened_by_phone == TRUE & wants_to_cancel == 1 ~ "Call: Cancel/Reschedule",
        intervened_by_phone == TRUE & wants_to_cancel == 0 ~ "Call: Chased & Confirmed"
      ),
      
      # Stage 3: Final Clinical Outcome
      Stage3_Outcome = final_attendance
    ) %>%
    # Count the volumes for each unique combination of paths
    count(Stage1_Arm, Stage2_Action, Stage3_Outcome, name = "Patients") %>%
    mutate(Total_Cohort = sum(Patients))
  
  # 2. Build the Alluvial Flow Chart
  ggplot(flow_data, aes(y = Patients, axis1 = Stage1_Arm, axis2 = Stage2_Action, axis3 = Stage3_Outcome)) +
    
    # Draw the flowing ribbons connecting the stages
    geom_alluvium(aes(fill = Stage3_Outcome), width = 1/4, alpha = 0.6, curve_type = "cubic") +
    
    # Draw the solid vertical segments (Strata)
    geom_stratum(width = 1/4, fill = "grey20", color = "white", alpha = 0.9) +
    
    # Add the text/stats labels inside the segments
    geom_text(stat = "stratum", 
              aes(label = after_stat(paste0(stratum, "\n(n=", count, ")"))), 
              color = "white", size = 3.5, fontface = "bold") +
    
    # Styling and scales
    scale_fill_manual(values = c("Attended" = "#2ecc71", "DNA" = "#e74c3c", "Cancelled" = "#95a5a6")) +
    scale_x_discrete(limits = c("1. Trial Arm", "2. Interactive Pathway", "3. Final Outcome"), expand = c(0.05, 0.05)) +
    labs(
      title = "Patient Flow & Attrition Pipeline",
      subtitle = "Tracing total patient volumes through automated text branches, coordinator queues, and final outcomes.",
      y = "Total Patient Volume",
      fill = "Final Attendance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      legend.position = "bottom"
    )
}


library(ggridges)

plot_branch_timing_ridge <- function(trial_data) {
  
  if (!"wants_to_cancel" %in% names(trial_data)) trial_data$wants_to_cancel <- 0
  
  timing_data <- trial_data %>%
    filter(final_attendance != "Censored") %>%
    mutate(
      branch_category = case_when(
        trial_arm == "Control" ~ "Control Arm",
        trial_arm == "Intervention" & intervened_by_phone == TRUE & wants_to_cancel == 1 ~ "Intervention: Cancel/Reschedule Call",
        trial_arm == "Intervention" & intervened_by_phone == TRUE & wants_to_cancel == 0 ~ "Intervention: Confirmed via Call",
        trial_arm == "Intervention" & intervened_by_phone == FALSE ~ "Intervention: Text / Timeout Track",
        TRUE ~ "Standard Low-Risk"
      ),
      # Calculate total hours from arrival to appointment
      total_journey_hours = (appt_time - arrival_time) / 60,
      # Time spent waiting/processing before action
      action_delay_hours = case_when(
        !is.na(call_time) ~ (call_time - arrival_time) / 60,
        TRUE ~ (tier2_text_time - arrival_time) / 60
      )
    )
  
  ggplot(timing_data, aes(x = total_journey_hours, y = branch_category, fill = branch_category)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
    scale_fill_viridis_d(option = "mako", alpha = 0.8) +
    labs(
      title = "Temporal Variation Across Pathway Branches",
      subtitle = "Distribution of total patient journey duration (Hours from List Extraction to Appointment) by branch",
      x = "Total Journey Duration (Hours)",
      y = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 10)
    )
}

plot_branch_milestone_spread <- function(trial_data, target_arm = "Intervention") {
  
  if (!"wants_to_cancel" %in% names(trial_data)) trial_data$wants_to_cancel <- 0
  
  step_data <- trial_data %>%
    filter(trial_arm == target_arm, final_attendance != "Censored") %>%
    mutate(
      pathway_type = case_when(
        intervened_by_phone == TRUE & wants_to_cancel == 1 ~ "Cancel / Reschedule Call",
        intervened_by_phone == TRUE & wants_to_cancel == 0 ~ "Confirmed via Call",
        intervened_by_phone == FALSE ~ "Text-Only / Timeout Path"
      ),
      # Compute cumulative hours from simulation start for each milestone
      t_arrival = arrival_time / 60,
      t_tier1   = tier1_text_time / 60,
      t_tier2   = tier2_text_time / 60,
      t_call    = ifelse(!is.na(call_time), call_time / 60, NA_real_),
      t_appt    = appt_time / 60
    ) %>%
    dplyr::select(patient_id, pathway_type, t_arrival, t_tier1, t_tier2, t_call, t_appt) %>%
    pivot_longer(cols = starts_with("t_"), names_to = "Milestone", values_to = "Hours") %>%
    filter(!is.na(Hours)) %>%
    mutate(
      Milestone = factor(Milestone, 
                         levels = c("t_arrival", "t_tier1", "t_tier2", "t_call", "t_appt"),
                         labels = c("1. List Extraction", "2. Tier 1 Text", "3. Tier 2 Touchpoint", "4. Coordinator Call", "5. Appointment"))
    )
  
  ggplot(step_data, aes(x = Milestone, y = Hours, fill = pathway_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("Milestone Timing Spread by Branch (", target_arm, " Arm)"),
      subtitle = "Shows median and variance (IQR) of absolute timeline hours for patients traveling down specific branches",
      x = "Pathway Milestone",
      y = "Simulation Time (Hours)",
      fill = "Branch Route"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
      legend.position = "bottom"
    )
}

plot_branch_latency_bars <- function(trial_data) {
  
  if (!"wants_to_cancel" %in% names(trial_data)) trial_data$wants_to_cancel <- 0
  
  latency_summary <- trial_data %>%
    filter(trial_arm == "Intervention", final_attendance != "Censored") %>%
    mutate(
      pathway_type = case_when(
        intervened_by_phone == TRUE & wants_to_cancel == 1 ~ "Cancel / Reschedule Call",
        intervened_by_phone == TRUE & wants_to_cancel == 0 ~ "Confirmed via Call",
        intervened_by_phone == FALSE ~ "Text-Only / Timeout Path"
      ),
      # Calculate specific operational lag intervals in hours
      lag_t1_to_t2 = (tier2_text_time - tier1_text_time) / 60,
      lag_t2_to_call = ifelse(!is.na(call_time), (call_time - tier2_text_time) / 60, NA_real_),
      lag_call_to_appt = ifelse(!is.na(call_time), (appt_time - call_time) / 60, (appt_time - tier2_text_time) / 60)
    ) %>%
    dplyr::select(pathway_type, lag_t1_to_t2, lag_t2_to_call, lag_call_to_appt) %>%
    pivot_longer(cols = starts_with("lag_"), names_to = "Transition", values_to = "Duration_Hours") %>%
    filter(!is.na(Duration_Hours) & Duration_Hours >= 0) %>%
    group_by(pathway_type, Transition) %>%
    summarise(
      mean_duration = mean(Duration_Hours),
      sd_duration = sd(Duration_Hours),
      .groups = "drop"
    ) %>%
    mutate(
      Transition_Label = case_when(
        Transition == "lag_t1_to_t2" ~ "Tier 1 Text -> Tier 2 Touchpoint",
        Transition == "lag_t2_to_call" ~ "Tier 2 Gate -> Coordinator Call",
        Transition == "lag_call_to_appt" ~ "Final Action -> Appointment"
      )
    )
  
  ggplot(latency_summary, aes(x = Transition_Label, y = mean_duration, fill = pathway_type)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_duration - sd_duration), ymax = mean_duration + sd_duration),
                  position = position_dodge(0.8), width = 0.25, color = "black") +
    scale_fill_brewer(palette = "Accent") +
    labs(
      title = "Operational Latency & Variance Between Pathway Steps",
      subtitle = "Mean duration (hours) and standard deviation error bars for each procedural transition",
      x = "Pathway Transition Step",
      y = "Mean Duration (Hours)",
      fill = "Branch Route"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"),
      legend.position = "bottom"
    )
}


plot_empirical_tracemap_static <- function(trial_data) {
  
  # Clean up the long pathway labels so they fit on the plot
  flow_data <- trial_data %>%
    mutate(
      path_clean = str_replace_all(pathway_outcome, "Intervention Arm: |Standard ", ""),
      path_clean = str_wrap(path_clean, width = 15)
    ) %>%
    count(risk_profile, trial_arm, path_clean, final_attendance, name = "Volume")
  
  ggplot(flow_data,
         aes(y = Volume, 
             axis1 = risk_profile, 
             axis2 = trial_arm, 
             axis3 = path_clean, 
             axis4 = final_attendance)) +
    geom_alluvium(aes(fill = final_attendance), width = 1/8, alpha = 0.7, color = "white") +
    geom_stratum(width = 1/4, fill = "grey90", color = "darkgray") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
    scale_x_discrete(limits = c("Risk Profile", "Trial Arm", "Pathway Outcome", "Final Attendance"), 
                     expand = c(0.1, 0.1)) +
    scale_fill_manual(values = c("Attended" = "#66c2a5", "DNA" = "#fc8d62", 
                                 "Cancelled" = "#8da0cb", "Censored" = "grey60")) +
    theme_minimal(base_size = 14) +
    labs(title = "Empirical Patient Flow Trace",
         subtitle = "Trace map generated directly from simulated volumes",
         y = "Number of Patients", x = "") +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}



generate_trial_process_suite <- function(trial_data) {
  
  # =========================================================================
  # 1. BUILD THE EVENT LOG
  # =========================================================================
  event_data <- trial_data %>%
    # Filter to actual trial participants to keep the visuals focused
    filter(trial_arm != "Not in Trial") %>% 
    dplyr::select(patient_id, trial_arm, risk_profile, sms_tier2_type, final_attendance,
           arrival_time, tier1_text_time, tier2_text_time, call_time, appt_time) %>%
    
    pivot_longer(
      cols = c(arrival_time, tier1_text_time, tier2_text_time, call_time, appt_time),
      names_to = "milestone",
      values_to = "time_mins",
      values_drop_na = TRUE
    ) %>%
    mutate(
      activity = case_when(
        milestone == "arrival_time"    ~ "1. Trial Entry",
        milestone == "tier1_text_time" ~ "2. Tier 1 Auto-SMS",
        milestone == "tier2_text_time" ~ paste0("3. Tier 2: ", sms_tier2_type),
        milestone == "call_time"       ~ "4. Coordinator Call",
        milestone == "appt_time"       ~ paste0("5. Outcome: ", final_attendance)
      ),
      # Convert simulation minutes to valid timestamps
      timestamp = as.POSIXct("2024-01-01 08:00:00") + (time_mins * 60),
      status = "complete",
      resource = "Clinic System",
      activity_instance = row_number()
    ) %>%
    arrange(patient_id, timestamp) %>%
    dplyr::select(patient_id, activity, timestamp, status, resource, activity_instance, trial_arm)
  
  trial_eventlog <- event_data %>%
    eventlog(
      case_id = "patient_id",
      activity_id = "activity",
      activity_instance_id = "activity_instance",
      lifecycle_id = "status",
      timestamp = "timestamp",
      resource_id = "resource"
    )
  
  # =========================================================================
  # 2. GENERATE THE CHARTS
  # =========================================================================
  
  # Chart 1: Frequency Trace Map (Patient Volumes)
  freq_map <- trial_eventlog %>% process_map(frequency("absolute"))
  
  # Chart 2: Performance Trace Map (Bottlenecks / Queue Delays)
  perf_map <- trial_eventlog %>% process_map(performance(mean, "hours"))
  
  # Chart 3: Trace Explorer (Visualizes every unique pathway permutation)
  # coverage = 1.0 means it will plot 100% of the variations in the data
  trace_exp <- trial_eventlog %>% trace_explorer(coverage = 1.0)
  
  # Chart 4: Animated Patient Flow
  # We color the tokens based on their trial arm so you can watch them diverge!
  animation_abs <- trial_eventlog %>%
    animate_process(
      mode = "absolute",
      mapping = token_aes(color = token_scale("trial_arm", 
                                              scale = "ordinal", 
                                              range = c("#3498db", "#e74c3c")))
    )
  animation_rel <- trial_eventlog %>%
    animate_process(
      mode = "relative",
      mapping = token_aes(color = token_scale("trial_arm", 
                                              scale = "ordinal", 
                                              range = c("#3498db", "#e74c3c")))
    )
  
  
  # =========================================================================
  # 3. RETURN AS NAMED LIST
  # =========================================================================
  return(list(
    event_log = trial_eventlog,
    freq_map = freq_map,
    perf_map = perf_map,
    trace_exp = trace_exp,
    animation_abs = animation_abs,
    animation_rel = animation_rel
  ))
}


# =========================================================================
# SECTION 4: MONTE CARLO PARAMETER STABILITY LOOP
# =========================================================================

run_stability_test <- function(iterations = 100, weeks_to_simulate = 26) {
  
  # Set your known ground-truth targets
  true_rr_tier2 <- 0.85 / 0.90 # Interactive (0.85) vs Passive Control (0.90)
  true_rr_tier3 <- (0.95 * 0.85 * 0.40) / (0.95 * 0.90) # Cumulative Tier 3 vs Control
  
  map_dfr(1:iterations, function(i) {
    # 1. Simulate a fresh trial dataset
    trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = weeks_to_simulate)
    
    # 2. Clean and build evaluation variables
    analysis_data <- trial_data %>%
      filter(!final_attendance %in% c("Censored", "Cancelled")) %>%
      mutate(
        dna_outcome = ifelse(final_attendance == "DNA", 1, 0),
        trial_arm = factor(trial_arm, levels = c("Control", "Intervention", "Not in Trial"))
      )
    
    high_risk_trial_cohorts <- analysis_data %>% 
      filter(trial_arm != "Not in Trial") %>%
      mutate(
        pp_exposure = case_when(
          trial_arm == "Control" ~ "1_Control", 
          trial_arm == "Intervention" & intervened_by_phone == FALSE & final_attendance != "Censored" ~ "2_Tier2_Interactive_Only",
          trial_arm == "Intervention" & intervened_by_phone == TRUE ~ "3_Tier3_Call_Reached",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(pp_exposure)) %>%
      mutate(
        pp_exposure = factor(pp_exposure, levels = c("1_Control", "2_Tier2_Interactive_Only", "3_Tier3_Call_Reached"))
      )
    
    # Safely fit the model (handles rare edge-case convergence warnings gracefully)
    model_fit <- tryCatch({
      glm(dna_outcome ~ pp_exposure + offset(log(ml_baseline_risk)),
          data = high_risk_trial_cohorts, family = poisson(link = "log"))
    }, error = function(e) { NULL })
    
    if (is.null(model_fit)) return(NULL)
    
    res <- coeftest(model_fit, vcov = sandwich)
    
    # Extract coefficients
    tibble(
      iteration = i,
      term = rownames(res),
      estimate = res[, "Estimate"],
      std_error = res[, "Std. Error"]
    )
  }) %>%
    mutate(
      Risk_Ratio = exp(estimate),
      Target_RR = case_when(
        grepl("Tier2", term) ~ true_rr_tier2,
        grepl("Tier3", term) ~ true_rr_tier3,
        TRUE ~ NA_real_
      )
    )
}

plot_stability_forest <- function(stability_results) {
  
  # 1. Calculate Summary Data for the Forest Plot
  summary_data <- stability_results %>%
    filter(!is.na(Target_RR)) %>%
    group_by(term) %>%
    summarise(
      Mean_Recovered_RR = mean(Risk_Ratio),
      # Calculate empirical 95% bounds based on the Monte Carlo distribution
      Conf_Low = quantile(Risk_Ratio, 0.025), 
      Conf_High = quantile(Risk_Ratio, 0.975), 
      Target_RR = first(Target_RR),
      .groups = "drop"
    )
  
  # 2. Build the Forest Plot
  p <- ggplot(summary_data, aes(y = term)) +
    
    # Ground-Truth Reference Lines
    geom_vline(aes(xintercept = Target_RR, color = term), 
               linetype = "dashed", linewidth = 1, alpha = 0.6) +
    
    # Recovered Point Estimates and 95% Simulation Intervals
    geom_pointrange(aes(x = Mean_Recovered_RR, 
                        xmin = Conf_Low, 
                        xmax = Conf_High, 
                        color = term),
                    size = 1.2, linewidth = 1.2) +
    
    # Safely constrain view without dropping wide Tier 3 intervals
    coord_cartesian(xlim = c(0.2, 1.2)) +
    
    labs(
      title = "Parameter stability: validation forest plot",
      subtitle = "Points represent mean recovered risk ratios with 95% simulation intervals.\nDashed lines indicate programmed ground-truth targets.",
      x = "Risk Ratio",
      y = "Exposure term"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold")
    )
  
  return(p)
}
plot_stability_forest(stability_results)

