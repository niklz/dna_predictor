# Load required libraries
library(tidyverse)
library(broom) # For tidy statistical models

set.seed(42) # For reproducibility

# 1. Generate a cohort of 1000 High-Risk Patients
n_patients <- 1000

evaluation_data <- tibble(
  appointment_id = 1:n_patients,
  risk_score = runif(n_patients, 0.75, 0.99), # All high risk
  # Randomly assign to Group A (Standard) or Group B (Cascade Policy)
  policy_arm = sample(c("Group A (Control)", "Group B (Cascade)"), n_patients, replace = TRUE)
)

# 2. Simulate the Intervention Funnel (Group B only)
evaluation_data <- evaluation_data %>%
  mutate(
    # TIER 2: 40% of Group B responds to the interactive text, 60% ignore it
    tier_2_status = case_when(
      policy_arm == "Group A (Control)" ~ "Not Attempted",
      TRUE ~ sample(c("Responded", "Ignored"), n(), prob = c(0.40, 0.60), replace = TRUE)
    ),
    
    # TIER 3: Only triggered if Group B ignored Tier 2
    tier_3_escalation = case_when(
      tier_2_status == "Ignored" ~ "Call Attempted",
      TRUE ~ "Not Attempted"
    ),
    
    # Did they answer the phone? (30% answer rate)
    tier_3_status = case_when(
      tier_3_escalation == "Call Attempted" ~ sample(c("Answered", "Voicemail"), n(), prob = c(0.30, 0.70), replace = TRUE),
      TRUE ~ "Not Attempted"
    )
  )

# 3. Simulate Ground Truth Outcomes based on interventions
# Operational Failure = "No-Show" or "Late Cancel"
# Operational Success = "Attended", "Timely Cancel", "Timely Reschedule"
evaluation_data <- evaluation_data %>%
  mutate(
    final_outcome = case_when(
      # Group A: High failure rate baseline (e.g., 40% failure)
      policy_arm == "Group A (Control)" ~ 
        sample(c("Success", "Failure"), n(), prob = c(0.60, 0.40), replace = TRUE),
      
      # Group B - Deflected at Tier 2: High success rate (e.g., 85% success)
      tier_2_status == "Responded" ~ 
        sample(c("Success", "Failure"), n(), prob = c(0.85, 0.15), replace = TRUE),
      
      # Group B - Reached at Tier 3: Good success rate (e.g., 75% success)
      tier_3_status == "Answered" ~ 
        sample(c("Success", "Failure"), n(), prob = c(0.75, 0.25), replace = TRUE),
      
      # Group B - Unreachable (Ignored text, didn't answer call): Very high failure rate
      tier_3_status == "Voicemail" ~ 
        sample(c("Success", "Failure"), n(), prob = c(0.50, 0.50), replace = TRUE)
    ),
    # Convert to binary numeric for modeling (1 = Failure/No-Show, 0 = Success)
    is_failure = if_else(final_outcome == "Failure", 1, 0)
  )


# Calculate the overall Failure Rate between the two arms
top_line_results <- evaluation_data %>%
  group_by(policy_arm) %>%
  summarise(
    total_patients = n(),
    failures = sum(is_failure),
    failure_rate = mean(is_failure)
  ) %>%
  mutate(failure_rate_pct = scales::percent(failure_rate, accuracy = 0.1))

print(top_line_results)

# Statistical Test: Logistic Regression to prove significance
# Are the odds of failure significantly lower in Group B?
adjusted_model <- glm(
  is_failure ~ policy_arm + risk_score, 
  data = evaluation_data, 
  family = "binomial"
)

# Tidy summary showing the true, adjusted impact of the Cascade Policy
library(broom)
tidy(adjusted_model, exponentiate = TRUE, conf.int = TRUE)



# Filter to only the patients who went through the ML workflow
cascade_data <- evaluation_data %>% filter(policy_arm == "Group B (Cascade)")

# 1. Tier 2 Deflection Rate (How much work did the automation save us?)
deflection_rate <- cascade_data %>%
  count(tier_2_status) %>%
  mutate(pct = scales::percent(n / sum(n)))

print("--- Tier 2 Automation Deflection ---")
print(deflection_rate)

# 2. Tier 3 Escalation Burden (How many calls do staff actually need to make?)
# Useful for capacity planning: "Out of N high-risk patients, expect to call Y of them."
escalation_burden <- cascade_data %>%
  count(tier_3_escalation) %>%
  mutate(pct = scales::percent(n / sum(n)))

print("--- Tier 3 Staff Burden ---")
print(escalation_burden)

# 3. The Unreachables (Where is the system failing?)
# Patients who ignored the text AND sent the call to voicemail
unreachable_rate <- cascade_data %>%
  filter(tier_3_status == "Voicemail") %>%
  summarise(
    unreachable_count = n(),
    pct_of_group_b = scales::percent(n() / nrow(cascade_data))
  )

print("--- The Unreachables ---")
print(unreachable_rate)
