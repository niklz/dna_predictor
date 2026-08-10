library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(broom)
library(sandwich)
library(lmtest)
library(tidyr)

# =========================================================================
# 1. CORE SIMULATION ENGINE
# =========================================================================
simulate_clinical_trial_advanced <- function(
    weeks_to_simulate = 8,
    total_patients_per_week = 1000,
    ml_high_risk_prop = 0.25,
    trial_allocation_ratio = 0.50,
    coordinator_capacity = 1,
    max_calls_per_week = 100,
    base_response_prob = 0.30, 
    base_prob_cancel = 0.15,
    rr_tier1          = 0.95,  
    rr_tier2_passive  = 0.90,  
    rr_tier2_interact = 0.85,  
    rr_tier3_call     = 0.40,
    cancellation_scalar = 1.25
) {
  
  mins_per_day <- 1440
  mins_per_week <- 7 * mins_per_day
  total_sim_time <- weeks_to_simulate * mins_per_week
  work_mins_per_week <- 5 * 8 * 60
  mean_handling_time <- work_mins_per_week / max_calls_per_week
  
  num_weeks <- weeks_to_simulate
  weekly_total_sizes <- rpois(num_weeks, lambda = total_patients_per_week)
  
  trial_manifest <- data.frame(week = rep(1:num_weeks, times = weekly_total_sizes)) %>%
    mutate(
      patient_id = paste0("pt_", row_number()),
      risk_profile = ifelse(runif(n()) < ml_high_risk_prop, "high_risk", "low_risk"),
      ml_baseline_risk = ifelse(risk_profile == "low_risk", 
                                rbeta(n(), shape1 = 2, shape2 = 50), 
                                rbeta(n(), shape1 = 5, shape2 = 20)),
      trial_arm = case_when(
        risk_profile == "low_risk" ~ "not_in_trial",
        risk_profile == "high_risk" & runif(n()) < trial_allocation_ratio ~ "control",
        TRUE ~ "intervention"
      )
    )
  
  intervention_df <- trial_manifest %>% filter(trial_arm == "intervention") %>% mutate(sim_idx = row_number())
  standard_df <- trial_manifest %>% filter(trial_arm != "intervention") %>% mutate(sim_idx = row_number())
  
  intervention_arrivals <- (intervention_df$week - 1) * mins_per_week
  standard_arrivals <- (standard_df$week - 1) * mins_per_week
  
  env <- simmer("trial_pathway")
  
  intervention_path <- trajectory("high_risk_intervention_journey") %>%
    set_attribute("arrival_time", function() { simmer::now(env) }) %>%
    set_attribute("is_intervention_arm", 1) %>%
    set_attribute("days_until_appt", function() { runif(1, min = 4, max = 14) }) %>%
    set_attribute("appt_time", function() { simmer::now(env) + (get_attribute(env, "days_until_appt") * mins_per_day) }) %>%
    set_attribute("tier1_text_time", function() { 
      base_t1 <- get_attribute(env, "appt_time") - (7 * mins_per_day) + rnorm(1, 0, 120)
      max(simmer::now(env), base_t1) 
    }) %>%
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
        
        patient_name <- simmer::get_name(env)
        current_idx <- as.numeric(sub(".*_", "", patient_name)) + 1
        patient_risk <- intervention_df$ml_baseline_risk[current_idx]
        
        dyn_response_prob <- base_response_prob * (1.5 * (1 - patient_risk))
        dyn_response_prob <- pmin(0.85, pmax(0.05, dyn_response_prob)) 
        
        did_respond <- rbinom(1, 1, dyn_response_prob)
        if (did_respond == 0 || response_t > cutoff_t) { return(3) } 
        
        is_intervention <- get_attribute(env, "is_intervention_arm") == 1 
        
        if (is_intervention) {
          dyn_cancel_prob <- base_prob_cancel * cancellation_scalar
        } else {
          dyn_cancel_prob <- base_prob_cancel
        }
        
        dyn_cancel_prob <- pmin(0.8, pmax(0.01, dyn_cancel_prob))
        return(sample(c(1, 2), size = 1, prob = c(1 - dyn_cancel_prob, dyn_cancel_prob)))
      },
      continue = c(TRUE, TRUE, TRUE),
      
      trajectory("intent_confirm") %>%
        set_attribute("response_type", 1) %>%
        set_attribute("outcome_status", 2),
      
      trajectory("intent_reschedule_cancel") %>%
        set_attribute("response_type", 2) %>%
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
      
      trajectory("intent_no_response") %>%
        set_attribute("response_type", 3) %>%
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
        set_attribute("wants_to_cancel", function() { 
          dyn_cancel_prob <- base_prob_cancel * cancellation_scalar
          dyn_cancel_prob <- pmin(0.8, pmax(0.01, dyn_cancel_prob))
          return(rbinom(1, 1, dyn_cancel_prob)) 
        }) %>%
        release("coordinator", 1) %>%
        set_attribute("outcome_status", 1)
    )
  
  standard_path <- trajectory("standard_care_journey") %>%
    set_attribute("arrival_time", function() { simmer::now(env) }) %>%
    set_attribute("is_intervention_arm", 0) %>%
    set_attribute("days_until_appt", function() { runif(1, min = 4, max = 14) }) %>%
    set_attribute("appt_time", function() { simmer::now(env) + (get_attribute(env, "days_until_appt") * mins_per_day) }) %>%
    set_attribute("tier1_text_time", function() { 
      base_t1 <- get_attribute(env, "appt_time") - (7 * mins_per_day) + rnorm(1, 0, 120)
      max(simmer::now(env), base_t1)
    }) %>%
    set_attribute("tier2_text_time", function() { 
      base_t2 <- get_attribute(env, "appt_time") - (3 * mins_per_day) + rnorm(1, 0, 60)
      max(get_attribute(env, "tier1_text_time"), base_t2) 
    })
  
  working_hours <- schedule(timetable = c(0, 9*60, 17*60), values = c(0, coordinator_capacity, 0), period = mins_per_day)
  
  env %>%
    add_resource("coordinator", capacity = working_hours) %>%
    add_generator("intervention_", intervention_path, at(intervention_arrivals), mon = 2) %>%
    add_generator("standard_", standard_path, at(standard_arrivals), mon = 2) %>%
    run(until = total_sim_time)
  
  sim_arrivals <- get_mon_arrivals(env) %>% 
    mutate(arm_group = sub("_.*", "", name), sim_idx = as.numeric(sub(".*_", "", name)) + 1)
  
  safe_max <- function(x) { if (all(is.na(x))) return(NA_real_) else return(max(x, na.rm = TRUE)) }
  
  sim_attributes <- get_mon_attributes(env) %>%
    group_by(name, key) %>% slice_tail(n = 1) %>% ungroup() %>%
    pivot_wider(names_from = key, values_from = value) %>%
    group_by(name) %>% summarise(across(everything(), safe_max), .groups = "drop")
  
  sim_processed <- sim_arrivals %>% left_join(sim_attributes, by = "name")
  
  expected_cols <- c("outcome_status", "call_time", "appt_time", "tier1_text_time", "tier2_text_time", "wants_to_cancel", "days_until_appt", "response_type")
  for (col in expected_cols) { if (!col %in% names(sim_processed)) sim_processed[[col]] <- NA_real_ }
  
  intervention_manifest <- intervention_df %>%
    left_join(sim_processed %>% filter(arm_group == "intervention"), by = "sim_idx") %>% 
    dplyr::select(-sim_idx, -name, -arm_group) %>% 
    mutate(outcome_status = replace_na(outcome_status, -1), sms_tier2_type = "interactive_gate")
  
  non_intervention_manifest <- standard_df %>%
    left_join(sim_processed %>% filter(arm_group == "standard"), by = "sim_idx") %>%
    dplyr::select(-sim_idx, -name, -arm_group) %>% 
    mutate(outcome_status = -2, wants_to_cancel = 0, sms_tier2_type = "passive_reminder", call_time = NA_real_)
  
  final_trial_dataset <- bind_rows(intervention_manifest, non_intervention_manifest) %>%
    mutate(
      intervened_by_phone = ifelse(outcome_status == 1, TRUE, FALSE),
      is_truncated = ifelse(appt_time > total_sim_time, TRUE, FALSE),
      pathway_outcome = case_when(
        is_truncated ~ "simulation_truncated",
        outcome_status == -1 ~ "intervention_pending_queue",
        trial_arm == "not_in_trial" ~ "standard_low_risk_track",
        trial_arm == "control" ~ "standard_high_risk_control",
        outcome_status == 0 ~ "intervention_queue_timeout",
        outcome_status == 1 ~ "intervention_processed_coordinator", 
        TRUE ~ "completed_auto_text_only"
      ),
      modulated_dna_prob = ml_baseline_risk * ifelse(appt_time > arrival_time, rr_tier1, 1) * ifelse(sms_tier2_type == "passive_reminder", rr_tier2_passive, rr_tier2_interact) * ifelse(pathway_outcome == "intervention_processed_coordinator", rr_tier3_call, 1),
      modulated_dna_prob = pmin(1, pmax(0, modulated_dna_prob)),
      modulated_cancel_prob = base_prob_cancel * ifelse(trial_arm == "intervention", cancellation_scalar, 1),
      modulated_cancel_prob = pmin(1, pmax(0, modulated_cancel_prob)),
      final_attendance = case_when(
        pathway_outcome %in% c("simulation_truncated", "intervention_pending_queue", "intervention_queue_timeout") ~ "censored",
        runif(n()) < modulated_cancel_prob ~ "cancelled",
        runif(n()) < modulated_dna_prob ~ "dna",
        TRUE ~ "attended"
      )
    )
  
  sop_dataset <- final_trial_dataset %>%
    dplyr::select(
      patient_id, trial_arm, risk_profile, ml_baseline_risk,        
      sms_tier2_type, intervened_by_phone, final_attendance,        
      arrival_time, appt_time, tier1_text_time, tier2_text_time, call_time, pathway_outcome, response_type 
    )
  return(sop_dataset)
}


# =========================================================================
# 2. PROCESS MINING & EXPLORATION VISUALIZATIONS
# =========================================================================

plot_trial_journeys <- function(trial_data, num_samples_per_group = 3) {
  
  sampled_patients <- trial_data %>%
    filter(!final_attendance %in% c("censored", "cancelled")) %>%
    group_by(trial_arm, final_attendance) %>%
    slice_sample(n = num_samples_per_group, replace = TRUE) %>%
    ungroup() %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    mutate(patient_plot_order = reorder(as.factor(patient_id), appt_time))
  
  touchpoint <- sampled_patients %>%
    transmute(
      patient_id = patient_plot_order, trial_arm,
      list_extraction = arrival_time / 60,
      tier1_automated_text = tier1_text_time / 60,
      tier2_passive_reminder = ifelse(sms_tier2_type == "passive_reminder", tier2_text_time / 60, NA_real_),
      tier2_interactive_gate = ifelse(sms_tier2_type == "interactive_gate", tier2_text_time / 60, NA_real_),
      coordinator_call = call_time / 60, 
      appointment = appt_time / 60
    ) %>%
    pivot_longer(cols = -c(patient_id, trial_arm), names_to = "milestone", values_to = "time_hours") %>%
    filter(!is.na(time_hours) & time_hours >= 0)
  
  lifespans <- sampled_patients %>%
    transmute(patient_id = patient_plot_order, trial_arm, final_attendance,
              start = pmax(0, pmin(arrival_time, tier1_text_time, na.rm = TRUE) / 60),
              end = appt_time / 60)
  
  max_hrs <- max(lifespans$end, na.rm = TRUE)
  weekend_starts <- seq(120, max_hrs, by = 168)
  weekends <- data.frame(xmin = weekend_starts, xmax = weekend_starts + 48) %>% filter(xmin <= max_hrs)
  
  ggplot() +
    geom_rect(data = weekends, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey90", alpha = 0.6, inherit.aes = FALSE) +
    geom_segment(data = lifespans, aes(x = start, xend = end, y = patient_id, yend = patient_id, color = final_attendance),
                 linewidth = 1.2, alpha = 0.75) +
    geom_point(data = touchpoint, aes(x = time_hours, y = patient_id, shape = milestone, fill = milestone),
               size = 3.2, stroke = 0.6, color = "black", alpha = 0.8) +
    scale_shape_manual(values = c("list_extraction"=21, "tier1_automated_text"=22, "tier2_passive_reminder"=23, 
                                  "tier2_interactive_gate"=24, "coordinator_call"=25, "appointment"=21)) +
    scale_fill_manual(values = c("list_extraction"="#66c2a5", "tier1_automated_text"="#fc8d62", 
                                 "tier2_passive_reminder"="#8da0cb", "tier2_interactive_gate"="#e78ac3", 
                                 "coordinator_call"="#a6d854", "appointment"="#ffd92f")) +
    scale_colour_manual(values = c("attended" = "darkolivegreen3", "dna" = "coral2")) +
    facet_wrap(~ trial_arm, scales = "free_y", ncol = 1) +
    labs(title = "Patient journey timelines", subtitle = "Grey bands indicate weekends.",
         x = "Time (hours from simulation start)", y = "") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", axis.text.y = element_blank())
}


analyze_pathway_volumes <- function(trial_data) {
  volume_summary <- trial_data %>%
    count(trial_arm, risk_profile, sms_tier2_type, intervened_by_phone, final_attendance) %>%
    mutate(percentage = round(n / sum(n) * 100, 2))
  
  print("--- Patient flow volume breakdown ---")
  print(volume_summary)
  return(volume_summary)
}


plot_transition_time_distributions <- function(trial_data) {
  transition_data <- trial_data %>%
    filter(!final_attendance %in% c("censored", "cancelled")) %>%
    mutate(
      dur_entry_to_t1 = (tier1_text_time - arrival_time) / 60,
      dur_t1_to_t2    = (tier2_text_time - tier1_text_time) / 60,
      dur_t2_to_call  = ifelse(intervened_by_phone == TRUE, (call_time - tier2_text_time) / 60, NA_real_),
      dur_call_to_appt = ifelse(intervened_by_phone == TRUE, (appt_time - call_time) / 60, (appt_time - tier2_text_time) / 60)
    ) %>%
    dplyr::select(patient_id, trial_arm, dur_entry_to_t1, dur_t1_to_t2, dur_t2_to_call, dur_call_to_appt) %>%
    pivot_longer(
      cols = starts_with("dur_"),
      names_to = "transition_stage",
      values_to = "hours"
    ) %>%
    filter(!is.na(hours) & hours >= 0) %>%
    mutate(
      stage_label = case_when(
        transition_stage == "dur_entry_to_t1"  ~ "1. Trial entry -> tier 1 text",
        transition_stage == "dur_t1_to_t2"     ~ "2. Tier 1 text -> tier 2 touchpoint",
        transition_stage == "dur_t2_to_call"   ~ "3. Tier 2 touchpoint -> coordinator call",
        transition_stage == "dur_call_to_appt" ~ "4. Final action -> appointment"
      )
    )
  
  ggplot(transition_data, aes(x = hours, fill = trial_arm)) +
    geom_density(alpha = 0.5, color = NA, adjust = 1.5) +
    facet_wrap(~ stage_label, scales = "free") +
    labs(
      title = "Distribution of transition times between pathway stages",
      subtitle = "Evaluates operational latency and variance across trial arms",
      x = "Duration (hours)",
      y = "Density",
      fill = "Trial arm"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.text = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 11)
    )
}


generate_trial_process_suite <- function(trial_data) {
  event_data <- trial_data %>%
    filter(trial_arm != "not_in_trial") %>% 
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
        milestone == "arrival_time"    ~ "1. Trial entry",
        milestone == "tier1_text_time" ~ "2. Tier 1 auto-SMS",
        milestone == "tier2_text_time" ~ paste0("3. Tier 2: ", tolower(sms_tier2_type)),
        milestone == "call_time"       ~ "4. Coordinator call",
        milestone == "appt_time"       ~ paste0("5. Outcome: ", tolower(final_attendance))
      ),
      timestamp = as.POSIXct("2024-01-01 08:00:00") + (time_mins * 60),
      status = "complete",
      resource = "Clinic System",
      activity_instance = row_number()
    ) %>%
    arrange(patient_id, timestamp) %>%
    dplyr::select(patient_id, activity, timestamp, status, resource, activity_instance, trial_arm)
  
  trial_eventlog <- event_data %>%
    eventlog(
      case_id = "patient_id", activity_id = "activity", activity_instance_id = "activity_instance",
      lifecycle_id = "status", timestamp = "timestamp", resource_id = "resource"
    )
  
  freq_map <- trial_eventlog %>% process_map(frequency("absolute"))
  perf_map <- trial_eventlog %>% process_map(performance(mean, "hours"))
  trace_exp <- trial_eventlog %>% trace_explorer(coverage = 1.0)
  
  animation_abs <- trial_eventlog %>%
    animate_process(mode = "absolute", mapping = token_aes(color = token_scale("trial_arm", scale = "ordinal", range = c("#3498db", "#e74c3c"))))
  
  animation_rel <- trial_eventlog %>%
    animate_process(mode = "relative", mapping = token_aes(color = token_scale("trial_arm", scale = "ordinal", range = c("#3498db", "#e74c3c"))))
  
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
# 3. STATISTICAL EVALUATION & VALIDATION 
# =========================================================================

plot_clean_dual_patchwork <- function(
    cancel_model_results,
    wasted_model_results,
    true_rr_cancel    = 1.25,           
    true_rr_tier2_dna = 0.85 / 0.90,  
    true_rr_tier3_dna = (0.95 * 0.85 * 0.40) / (0.95 * 0.90)
) {
  
  cancel_plot_data <- tidy(cancel_model_results) %>%
    filter(term == "trial_armintervention") %>%
    mutate(
      risk_ratio = exp(estimate),
      conf_low   = exp(estimate - 1.96 * std.error),
      conf_high  = exp(estimate + 1.96 * std.error),
      tier = "Intervention arm effect",
      true_parameter = true_rr_cancel,
      term_label = paste0("Intervention vs control\n(true rr = ", round(true_rr_cancel, 3), ")")
    )
  
  p1 <- ggplot(cancel_plot_data, aes(x = risk_ratio, y = term_label)) +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), color = "#2c3e50", size = 0.8, linewidth = 1.2) +
    geom_point(aes(x = true_parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    labs(title = "Intervention impact on advanced cancellation", x = "Recovered risk ratio", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  wasted_plot_data <- tidy(wasted_model_results) %>%
    filter(term %in% c("pp_exposuretier2_interactive_only", "pp_exposuretier3_call_reached")) %>%
    mutate(
      risk_ratio = exp(estimate),
      conf_low   = exp(estimate - 1.96 * std.error),
      conf_high  = exp(estimate + 1.96 * std.error),
      tier = case_when(
        grepl("tier2", term) ~ "Tier 2 impacts",
        grepl("tier3", term) ~ "Tier 3 impacts"
      ),
      true_parameter = case_when(
        grepl("tier2", term) ~ true_rr_tier2_dna,
        grepl("tier3", term) ~ true_rr_tier3_dna
      ),
      term_label = case_when(
        grepl("tier2", term) ~ paste0("Interactive text vs passive text\n(true rr = ", round(true_rr_tier2_dna, 3), ")"),
        grepl("tier3", term) ~ paste0("Call vs text only\n(true rr = ", round(true_rr_tier3_dna, 3), ")")
      )
    )
  
  p2 <- ggplot(wasted_plot_data, aes(x = risk_ratio, y = term_label)) +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), color = "#2c3e50", size = 0.8, linewidth = 1.2) +
    geom_point(aes(x = true_parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    facet_grid(tier ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Intervention impact on DNA", x = "Recovered risk ratio", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  return((p1 / p2) + plot_layout(heights = c(1/3, 2/3)) + plot_annotation(title = "Validation of parameter recovery"))
}


run_stability <- function(iterations = 100, weeks_to_simulate = 26) {
  
  # static programmed parameters (per-protocol)
  true_rr_cancel <- 1.25
  true_rr_tier2 <- 0.85 / 0.90 
  true_rr_tier3 <- (0.95 * 0.85 * 0.40) / (0.95 * 0.90) 
  
  map_dfr(1:iterations, function(i) {
    # 1. simulate
    trial_data <- simulate_clinical_trial_advanced(weeks_to_simulate = weeks_to_simulate)
    
    # 2. clean and build evaluation variables
    analysis_data <- trial_data %>%
      filter(final_attendance != "censored") %>%
      mutate(
        wasted_slot_outcome = ifelse(final_attendance == "dna", 1, 0),
        cancellation_outcome = ifelse(final_attendance == "cancelled", 1, 0),
        trial_arm = factor(trial_arm, levels = c("control", "intervention", "not_in_trial")),
        pp_exposure = case_when(
          trial_arm == "control" ~ "control", 
          trial_arm == "intervention" & intervened_by_phone == FALSE ~ "tier2_interactive_only",
          trial_arm == "intervention" & intervened_by_phone == TRUE ~ "tier3_call_reached",
          TRUE ~ NA_character_
        ),
        pp_exposure = factor(pp_exposure, levels = c("control", "tier2_interactive_only", "tier3_call_reached"))
      )
    
    # 3. DYNAMIC ITT EXPECTATION CALCULATION
    # empirical proportions for the intervention arm
    int_prop <- analysis_data %>%
      filter(trial_arm == "intervention") %>%
      summarise(
        p_cancel = mean(final_attendance == "cancelled"),
        p_tier2  = mean(final_attendance != "cancelled" & intervened_by_phone == FALSE),
        p_tier3  = mean(final_attendance != "cancelled" & intervened_by_phone == TRUE)
      )
    
    # empirical proportions for the control arm
    ctrl_prop <- analysis_data %>%
      filter(trial_arm == "control") %>%
      summarise(
        p_cancel = mean(final_attendance == "cancelled"),
        p_active = mean(final_attendance != "cancelled")
      )
    
    # calculate the emergent weighted average modifier for this specific run
    eff_numerator   <- (int_prop$p_cancel * 0) + (int_prop$p_tier2 * 0.85) + (int_prop$p_tier3 * (0.85 * 0.40))
    eff_denominator <- (ctrl_prop$p_cancel * 0) + (ctrl_prop$p_active * 0.90)
    
    dynamic_itt_dna <- eff_numerator / eff_denominator
    
    # 4. run models
    cancel_data <- analysis_data %>% filter(trial_arm %in% c("control", "intervention"))
    cancel_fit <- tryCatch({ glm(cancellation_outcome ~ trial_arm + ml_baseline_risk, data = cancel_data, family = poisson(link = "log")) }, error = function(e) { NULL })
    cancel_res <- if (!is.null(cancel_fit)) coeftest(cancel_fit, vcov = sandwich) else NULL
    
    global_dna_data <- analysis_data %>% filter(trial_arm %in% c("control", "intervention"))
    global_fit <- tryCatch({ glm(wasted_slot_outcome ~ trial_arm + ml_baseline_risk, data = global_dna_data, family = poisson(link = "log")) }, error = function(e) { NULL })
    global_res <- if (!is.null(global_fit)) coeftest(global_fit, vcov = sandwich) else NULL
    
    wasted_data <- analysis_data %>% filter(trial_arm %in% c("control", "intervention") & final_attendance != "cancelled")
    wasted_fit <- tryCatch({ glm(wasted_slot_outcome ~ pp_exposure + ml_baseline_risk, data = wasted_data, family = poisson(link = "log")) }, error = function(e) { NULL })
    wasted_res <- if (!is.null(wasted_fit)) coeftest(wasted_fit, vcov = sandwich) else NULL
    
    # 5. bind rows and inject targets directly
    bind_rows(
      if (!is.null(cancel_res)) tibble(
        iteration = i, term = rownames(cancel_res), 
        estimate = cancel_res[, "Estimate"], std_error = cancel_res[, "Std. Error"], 
        endpoint = "cancellation", 
        target_rr = ifelse(term == "trial_armintervention", true_rr_cancel, NA_real_)
      ) else NULL,
      
      if (!is.null(global_res)) tibble(
        iteration = i, term = rownames(global_res), 
        estimate = global_res[, "Estimate"], std_error = global_res[, "Std. Error"], 
        endpoint = "global_itt_dna", 
        target_rr = ifelse(term == "trial_armintervention", dynamic_itt_dna, NA_real_)
      ) else NULL,
      
      if (!is.null(wasted_res)) tibble(
        iteration = i, term = rownames(wasted_res), 
        estimate = wasted_res[, "Estimate"], std_error = wasted_res[, "Std. Error"], 
        endpoint = "pp_dna",
        target_rr = case_when(
          grepl("tier2", term) ~ true_rr_tier2,
          grepl("tier3", term) ~ true_rr_tier3,
          TRUE ~ NA_real_
        )
      ) else NULL
    )
  }) %>%
    mutate(risk_ratio = exp(estimate))
}

plot_stability <- function(stability_results) {
  
  # 1. summarize data and extract the targets automatically
  summary_data <- stability_results %>%
    group_by(endpoint, term) %>%
    summarise(
      risk_ratio = exp(mean(estimate, na.rm = TRUE)),
      conf_low   = exp(mean(estimate - 1.96 * std_error, na.rm = TRUE)), 
      conf_high  = exp(mean(estimate + 1.96 * std_error, na.rm = TRUE)), 
      true_parameter = mean(target_rr, na.rm = TRUE), # dynamic extraction
      .groups    = "drop"
    )
  
  # panel 1: advance cancellations (ITT)
  cancel_plot_data <- summary_data %>%
    filter(endpoint == "cancellation" & grepl("intervention", tolower(term))) %>%
    mutate(term_label = paste0("Intervention vs control\n(true RR = ", round(true_parameter, 3), ")"))
  
  p1 <- ggplot(cancel_plot_data, aes(x = risk_ratio, y = term_label)) +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), color = "#2c3e50", size = 0.8, linewidth = 1.2) +
    geom_point(aes(x = true_parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    labs(title = "Macro impact: advance cancellations (ITT)", x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  # panel 2: global DNA (ITT)
  global_dna_data <- summary_data %>%
    filter(endpoint == "global_itt_dna" & grepl("intervention", tolower(term))) %>%
    mutate(term_label = paste0("Intervention vs control\n(true RR = ", round(true_parameter, 3), ")"))
  
  p2 <- ggplot(global_dna_data, aes(x = risk_ratio, y = term_label)) +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), color = "#2c3e50", size = 0.8, linewidth = 1.2) +
    geom_point(aes(x = true_parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    coord_cartesian(xlim = c(0.2, 1.2)) + 
    labs(title = "Global impact: total DNAs (ITT)", x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  # panel 3: mechanism (PP)
  wasted_plot_data <- summary_data %>%
    filter(endpoint == "pp_dna" & grepl("tier2|tier3", tolower(term))) %>%
    mutate(
      tier = case_when(
        grepl("tier2", tolower(term)) ~ "Tier 2 impacts",
        grepl("tier3", tolower(term)) ~ "Tier 3 impacts"
      ),
      term_label = case_when(
        grepl("tier2", tolower(term)) ~ paste0("Interactive text vs passive text\n(true RR = ", round(true_parameter, 3), ")"),
        grepl("tier3", tolower(term)) ~ paste0("Call vs text only\n(true RR = ", round(true_parameter, 3), ")")
      )
    )
  
  p3 <- ggplot(wasted_plot_data, aes(x = risk_ratio, y = term_label)) +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), color = "#2c3e50", size = 0.8, linewidth = 1.2) +
    geom_point(aes(x = true_parameter), color = "#e74c3c", size = 4, shape = 18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
    coord_cartesian(xlim = c(0.2, 1.2)) + 
    facet_grid(tier ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Operational mechanism: DNAs by exposure (PP)", x = "Recovered risk ratio", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face="bold"))
  
  # combine using patchwork
  final_plot <- (p2 / p1 / p3) + 
    plot_layout(heights = c(1, 1, 2)) + 
    plot_annotation(
      title = "Evaluation model stability over 100 replications"
    )
  
  return(final_plot)
}
