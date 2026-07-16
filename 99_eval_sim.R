library(simmer)
library(simmer.plot)
library(dplyr)
library(ggplot2)
library(tidyr)

simulate_pathway <- function(
    weeks_to_simulate = 8,
    patients_per_week = 250,
    coordinator_capacity = 1,      # 1 FTE
    max_calls_per_week = 100,      # Max real-world throughput per FTE
    response_prob = 0.30,          # 30% reply to interactive text
    prob_wish_to_cancel = 0.20,    # Of those reached, 20% want to cancel/rearrange
    baseline_dna_rate = 0.15,      # 15% DNA if not intervened
    confirmed_dna_rate = 0.05      # 5% DNA if they verbally confirmed
) {
  
  # Time definitions (1 unit = 1 minute)
  mins_per_day <- 1440
  mins_per_week <- 7 * mins_per_day
  total_sim_time <- weeks_to_simulate * mins_per_week
  interarrival_time <- mins_per_week / patients_per_week
  
  # Calculate realistic service time per patient task to cap throughput.
  # 1 FTE has 2400 working minutes/week. To cap at 100 calls, 
  # each patient task must take 24 minutes of coordinator resource time.
  work_mins_per_week <- 5 * 8 * 60 # 5 days * 8 hours = 2400 mins
  mean_handling_time <- work_mins_per_week / max_calls_per_week
  
  # Initialize simulation environment
  env <- simmer("Outpatient_Pathway")
  
  # 1. Define the Patient Trajectory
  patient_path <- trajectory("Patient_Journey") %>%
    
    # Assign baseline attributes at T-7 days
    set_attribute("arrival_time", function() { now(env) }) %>%
    set_attribute("appt_time", function() { now(env) + (7 * mins_per_day) }) %>%
    set_attribute("responds_to_text", function() { rbinom(1, 1, response_prob) }) %>%
    
    # Advance time to T-3 days (Interactive Text sent)
    timeout(4 * mins_per_day) %>%
    
    # Reneging: If appt time passes while waiting in queue, they drop out of queue
    renege_in(
      t = function() { max(0, get_attribute(env, "appt_time") - now(env)) },
      out = trajectory() %>% set_attribute("outcome_status", 0) # 0 = Missed/Expired in Queue
    ) %>%
    
    # Enter the coordinator queue
    seize("coordinator", 1) %>%
    renege_abort() %>% # Safely reached before appt passed
    
    # Record when the coordinator call actually happened
    set_attribute("call_time", function() { now(env) }) %>%
    
    # Coordinator processes the patient (talk time + admin/logging time)
    timeout(function() { max(5, rnorm(1, mean = mean_handling_time, sd = 5)) }) %>%
    
    # Determine the patient's intent during the call
    set_attribute("wants_to_cancel", function() { rbinom(1, 1, prob_wish_to_cancel) }) %>%
    
    # Release coordinator and log successful contact
    release("coordinator", 1) %>%
    set_attribute("outcome_status", 1) # 1 = Reached by coordinator
  
  # 2. Define the Coordinator Schedule (Mon-Fri, 09:00 - 17:00)
  working_hours <- schedule(
    timetable = c(0, 9*60, 17*60), 
    values = c(0, coordinator_capacity, 0), 
    period = mins_per_day
  )
  
  plot(patient_path)
  
  # 3. Build and Run the Environment
  env %>%
    add_resource("coordinator", capacity = working_hours) %>%
    add_generator("Patient", patient_path, function() { rexp(1, 1/interarrival_time) }, mon = 2) %>%
    run(until = total_sim_time)
  
  # 4. Extract and Process the Data
  arrivals_data <- get_mon_arrivals(env)
  attributes_data <- get_mon_attributes(env) %>%
    group_by(name, key) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(name, key, value) %>%
    pivot_wider(names_from = key, values_from = value)
  
  # Ensure columns exist
  cols_to_check <- c("outcome_status", "call_time", "appt_time", "arrival_time", "responds_to_text", "wants_to_cancel")
  for (col in cols_to_check) {
    if (!col %in% names(attributes_data)) attributes_data[[col]] <- NA_real_
  }
  
  # Join and calculate final endpoints
  final_data <- arrivals_data %>%
    left_join(attributes_data, by = "name") %>%
    mutate(
      outcome_status = replace_na(as.numeric(outcome_status), -1),
      wants_to_cancel = replace_na(as.numeric(wants_to_cancel), 0),
      reached_in_time = ifelse(outcome_status == 1, TRUE, FALSE),
      hours_notice = (as.numeric(appt_time) - as.numeric(call_time)) / 60,
      
      # Categorize Pathway Outcome
      pathway_outcome = case_when(
        outcome_status == -1 ~ "Pending in Simulation",
        !reached_in_time ~ "Missed (Appt Passed in Queue)",
        reached_in_time & wants_to_cancel == 1 & hours_notice >= 48 ~ "Timely Cancel/Rearrange",
        reached_in_time & wants_to_cancel == 1 & hours_notice < 48 ~ "Late Cancel/Rearrange",
        reached_in_time & wants_to_cancel == 0 ~ "Verbally Confirmed Attendance"
      ),
      
      # Determine Actual Attendance
      final_attendance = case_when(
        pathway_outcome == "Pending in Simulation" ~ "Pending in Simulation",
        pathway_outcome == "Timely Cancel/Rearrange" ~ "Cancelled/Rearranged",
        pathway_outcome == "Late Cancel/Rearrange" ~ "Cancelled/Rearranged",
        
        # If they verbally confirmed, they have a lower DNA rate
        pathway_outcome == "Verbally Confirmed Attendance" & runif(n()) < confirmed_dna_rate ~ "DNA",
        pathway_outcome == "Verbally Confirmed Attendance" ~ "Attended",
        
        # If they were missed, they retain the baseline DNA rate
        pathway_outcome == "Missed (Appt Passed in Queue)" & runif(n()) < baseline_dna_rate ~ "DNA",
        TRUE ~ "Attended"
      )
    )
  
  return(final_data)
}

# Run the simulation for 8 weeks
sim_results <- simulate_pathway(weeks_to_simulate = 8)

# 1. See what happened on the pathway
table(sim_results$pathway_outcome)

# 2. See the final attendance states (including DNAs)
table(sim_results$final_attendance)


plot_patient_journeys <- function(sim_data, num_samples_per_group = 3) {
  
  # 1. Sample and clean the data
  sampled_patients <- sim_data %>%
    filter(pathway_outcome != "Pending in Simulation") %>%
    group_by(pathway_outcome) %>%
    slice_sample(n = num_samples_per_group) %>%
    ungroup() %>%
    # CRITICAL FIX 1: Convert to factor and drop unused levels
    mutate(name = droplevels(factor(name))) %>%
    # CRITICAL FIX 2: Reorder patient IDs chronologically by arrival time
    mutate(patient_id = reorder(name, arrival_time))
  
  # 2. Build milestone dataset
  milestones <- sampled_patients %>%
    transmute(
      patient_id = patient_id,
      outcome = pathway_outcome,
      Arrival = arrival_time / 60,
      Text_Sent = (arrival_time + (4 * 1440)) / 60,
      Call_Attempt = call_time / 60,
      Appointment = appt_time / 60
    ) %>%
    pivot_longer(
      cols = c(Arrival, Text_Sent, Call_Attempt, Appointment),
      names_to = "Milestone",
      values_to = "Time_Hours"
    ) %>%
    filter(!is.na(Time_Hours))
  
  # 3. Create lifespan segments
  lifespans <- sampled_patients %>%
    transmute(
      patient_id = patient_id,
      outcome = pathway_outcome,
      start = arrival_time / 60,
      end = appt_time / 60
    )
  
  # 4. Plot
  ggplot() +
    geom_segment(
      data = lifespans,
      aes(x = start, xend = end, y = patient_id, yend = patient_id, color = outcome),
      linewidth = 1.2, alpha = 0.4
    ) +
    geom_point(
      data = milestones,
      aes(x = Time_Hours, y = patient_id, shape = Milestone, fill = Milestone),
      size = 3.5, stroke = 1, color = "black"
    ) +
    scale_shape_manual(values = c("Arrival" = 21, "Text_Sent" = 23, "Call_Attempt" = 22, "Appointment" = 24)) +
    scale_fill_manual(values = c("Arrival" = "#999999", "Text_Sent" = "#56B4E9", "Call_Attempt" = "#E69F00", "Appointment" = "#D55E00")) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Visualizing Individual Patient Journeys",
      subtitle = "Ordered chronologically by arrival time",
      x = "Simulation Time (Hours)",
      y = "Individual Patients (Chronological order)",
      color = "Final Pathway Outcome",
      shape = "Key Event Milestones",
      fill = "Key Event Milestones"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = element_blank(),
      # CRITICAL FIX 3: Hide y-axis text to remove name clutter
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}


plot_patient_journeys(sim_results, num_samples_per_group = 100)
