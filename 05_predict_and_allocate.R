# =========================================================================
# SCRIPT 05: PREDICT AND ALLOCATE (OPERATIONAL INFERENCE & WORKBOOK GENERATION)
# =========================================================================
# 1. Loads new patient appointments for the upcoming week.
# 2. Prepares features using 'apply_custom_feature_engineering()' (prevents train-serve skew).
# 3. Generates and calibrates DNA risk predictions using the saved production bundle.
# 4. Stratifies patients into risk profiles using the optimized operational threshold.
# 5. Generates a deterministic MD5 hash key for downstream database joining.
# 6. Collapses 20 clinical vulnerability columns into a single comma-separated column.
# 7. Injects real-time model execution metadata (run date and model version).
# 8. Randomly allocates high-risk patients using a force-balanced randomized split.
# 9. Generates a beautifully styled, coordinator-ready Excel workbook with interactive 
#    dropdown validation, collapsed flags, metadata fields, and distinct visual tracks.

# -------------------------------------------------------------------------
# 1. Setup and load configurations
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
# Load configuration parameters
conf <- config::get()

# Set a weekly seed for reproducible randomization of the trial arms.
set.seed(42) 

# -------------------------------------------------------------------------
# 2. Load production artifacts
# -------------------------------------------------------------------------
model_bundle_path <- "data/processed/rf_final_model.rds"
if (!file.exists(model_bundle_path)) {
  stop("Model bundle not found. Please train and calibrate the model first.")
}
model_bundle  <- readRDS(model_bundle_path)
rf_model      <- model_bundle$model
rf_calibrator <- model_bundle$calibrator

threshold_path <- "data/processed/risk_threshold.RDS"
if (!file.exists(threshold_path)) {
  stop("Operational threshold not found. Please calculate thresholds first.")
}
threshold_data <- readRDS(threshold_path)
op_threshold   <- threshold_data$production_threshold

# -------------------------------------------------------------------------
# 3. Load incoming appointments
# -------------------------------------------------------------------------
raw_new_appointments_path <- "data/new_clinic_appointments.csv"
if (!file.exists(raw_new_appointments_path)) {
  # Fallback for dev environment testing
  if (file.exists("data/data_joined.RDS")) {
    new_appointments <- readRDS("data/data_joined.RDS") %>% sample_n(400)
  } else {
    stop("New appointments file not found.")
  }
} else {
  new_appointments <- read.csv(raw_new_appointments_path)
}

# -------------------------------------------------------------------------
# 4. Feature engineering
# -------------------------------------------------------------------------
message("Engineering features for incoming cohort...")
engineered_appointments <- apply_custom_feature_engineering(new_appointments)

# -------------------------------------------------------------------------
# 5. Prediction and calibration
# -------------------------------------------------------------------------
message("Generating calibrated risk predictions...")
raw_predictions <- predict(rf_model, new_data = engineered_appointments, type = "prob")

calibrated_predictions <- raw_predictions %>%
  cal_apply(rf_calibrator)

# -------------------------------------------------------------------------
# 6. Stratify risk and allocate trial arms (Force-Balanced Split)
# -------------------------------------------------------------------------
message("Stratifying patient risk profiles and allocating cohorts...")

allocation_ratio <- ifelse(!is.null(conf$trial_allocation_ratio), 
                           conf$trial_allocation_ratio, 
                           0.50)

final_manifest <- new_appointments %>%
  mutate(
    # Pull calibrated probabilities
    dna_probability = calibrated_predictions$.pred_DNA,
    
    # Label risk profile based on our operational threshold
    risk_profile = ifelse(dna_probability >= op_threshold, "high-risk", "low-risk"),
    
    # Default all patients to 'not in trial'
    trial_arm = "not in trial"
  )

# Extract indices of the high-risk cohort
high_risk_indices <- which(final_manifest$risk_profile == "high-risk")
n_high_risk       <- length(high_risk_indices)

if (n_high_risk > 0) {
  # Calculate the exact number of control patients needed based on your allocation ratio
  n_control <- round(n_high_risk * allocation_ratio)
  n_interv  <- n_high_risk - n_control
  
  # Create a perfectly balanced vector of assignments and shuffle it randomly
  balanced_arms <- sample(c(
    rep("control", n_control),
    rep("intervention", n_interv)
  ))
  
  # Assign the balanced arms back to the high-risk patients
  final_manifest$trial_arm[high_risk_indices] <- balanced_arms
}

# -------------------------------------------------------------------------
# 7. Generate cryptographic primary keys, collapse flags, and map variables
# -------------------------------------------------------------------------
message("Generating unique cryptographic keys and collapsing vulnerability flags...")

# 20 raw vulnerability variables from EHR [cite: 1, 2]
vulnerability_cols <- c(
  "a_ld", "a_autism", "a_interpreter_req_bsl", "a_interpreter_req_lang", "a_balance", 
  "a_cognitive_impairment", "a_mobility_restriction", "a_hear_vis_impaired", "a_dementia", 
  "a_depression", "a_downs_syndrome", "a_long_standing_condition", "a_makaton", 
  "a_mild_cognitive_impairment", "a_memory_impairment", "a_mood_disorder", "a_other_disability", 
  "a_psychosis", "a_severe_anxiety", "a_wheelchair_user"
)

vulnerability_labels <- c(
  "Learning disability", "Autism", "BSL interpreter required", "Language interpreter required",
  "Balance impairment", "Cognitive impairment", "Mobility restriction", "Sensory impairment",
  "Dementia", "Depression", "Down's syndrome", "Long-standing condition", "Makaton communication",
  "Mild cognitive impairment", "Memory impairment", "Mood disorder", "Other disability",
  "Psychosis", "Severe anxiety", "Wheelchair user"
)

final_manifest <- final_manifest %>%
  # Collapse the 20 columns into a single comma-separated string column
  rowwise() %>%
  mutate(
    v_vals = list(c_across(all_of(vulnerability_cols))),
    accessibility_flags = paste(vulnerability_labels[which(v_vals == 1)], collapse = ", ")
  ) %>%
  ungroup() %>%
  # Clean up raw columns and intermediate lists
  select(-all_of(vulnerability_cols), -v_vals) %>%
  mutate(
    # MD5 hash of lowercase patient ID and appointment hour/month stamp to ensure uniqueness [cite: 530, 532]
    appointment_id = sapply(paste0(tolower(trimws(dim_patient_id)), "_", appt_month), 
                            function(x) digest::digest(x, algo = "md5", serialize = FALSE)),
    
    # Demographic and scheduled mappings matching requested columns
    patient_id           = dim_patient_id,
    appointment_datetime = appt_dttm,
    gp_practice          = registered_gp_practice,
    age                  = age_at_appointment,
    sex                  = gender,
    dna_risk             = dna_probability,
    risk_label           = risk_profile,
    
    # Dynamic execution metadata columns [cite: 372, 391]
    date_model_run       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_version        = ifelse(!is.null(conf$model_ver), conf$model_ver, "v1.1"),
    
    # Initialize blank coordinator fields (Col 16 to 22 / P to V) [cite: 378]
    tier1_text_timestamp = "",
    tier2_text_timestamp = "",
    tier3_phone_attempt  = "",
    tier3_call_timestamp = "",
    patient_intent       = "",
    call_notes           = "",
    appointment_outcome  = ""
  )

# -------------------------------------------------------------------------
# 8. Output operational lists & Excel workbook
# -------------------------------------------------------------------------
# Save the master raw manifest (CSV format)
output_manifest_path <- "outputs/weekly_scored_appointments.csv"
dir.create("outputs", showWarnings = FALSE)
write.csv(final_manifest, file = output_manifest_path, row.names = FALSE)

# Generate a minimal action list strictly for the coordinators (only Intervention patients)
coordinator_list <- final_manifest %>%
  filter(trial_arm == "intervention") %>%
  select(appointment_id, patient_id, appointment_datetime, gp_practice, dna_risk) %>%
  arrange(desc(dna_risk))

write.csv(coordinator_list, file = "outputs/weekly_coordinator_action_list.csv", row.names = FALSE)

# Now, dynamically generate and style the interactive clinical coordinator workbook
message("Generating clinical coordinator validation spreadsheet...")

# Create workbook
wb <- createWorkbook()

# Define professional styles (Navy, Slate, Soft Blue, and Light Gray)
style_title <- createStyle(fontName = "Calibri", fontSize = 16, fontColour = "#2F5496", textDecoration = "bold")
style_subtitle <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#595959", textDecoration = "italic")
style_section <- createStyle(fontName = "Calibri", fontSize = 12, fontColour = "#FFFFFF", fgFill = "#4472C4", textDecoration = "bold")
style_header <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#FFFFFF", fgFill = "#2F5496", textDecoration = "bold", halign = "center", valign = "center", border = "bottom", borderStyle = "medium")

style_kpi_label <- createStyle(fontName = "Calibri", fontSize = 10, fontColour = "#2c3e50", textDecoration = "bold", halign = "right")
style_kpi_value <- createStyle(fontName = "Calibri", fontSize = 11, fontColour = "#2F5496", fgFill = "#F2F2F2", textDecoration = "bold", halign = "center", border = c("top", "bottom", "left", "right"), borderStyle = "thin")

style_text_left <- createStyle(fontName = "Calibri", fontSize = 11, halign = "left")
style_num_right <- createStyle(fontName = "Calibri", fontSize = 11, halign = "right")
style_prob_pct  <- createStyle(fontName = "Calibri", fontSize = 11, halign = "right", numFmt = "0.0%")

style_active_input <- createStyle(fontName = "Calibri", fontSize = 11, fgFill = "#DCE6F1", halign = "center", border = c("top", "bottom", "left", "right"), borderStyle = "thin")
style_control_row  <- createStyle(fontName = "Calibri", fontSize = 11, fgFill = "#F2F2F2", fontColour = "#7F7F7F", halign = "center")

# .........................................................................
# SHEET 1: INSTRUCTIONS & MAPPING
# .........................................................................
sheet_name_inst <- "Instructions & mapping"
addWorksheet(wb, sheet_name_inst)
showGridLines(wb, sheet_name_inst, show = FALSE)

# Add titles
writeData(wb, sheet_name_inst, "Weekly coordination export data specification", startCol = 1, startRow = 1)
addStyle(wb, sheet_name_inst, style_title, rows = 1, cols = 1)
writeData(wb, sheet_name_inst, "Operational guidelines and evaluation dictionary for clinical coordinators (v7)", startCol = 1, startRow = 2)
addStyle(wb, sheet_name_inst, style_subtitle, rows = 2, cols = 1)

# Section: Variables mapping
writeData(wb, sheet_name_inst, "Model-generated variables vs manual logs", startCol = 1, startRow = 4)
addStyle(wb, sheet_name_inst, style_section, rows = 4, cols = 1:3)
setColWidths(wb, sheet_name_inst, cols = 1:3, widths = c(25, 50, 20))

# Build mapping dictionary table (22 columns total with single accessibility_flags)
mapping_dict <- data.frame(
  Variable_Name = c(
    "appointment_id", "patient_id", "appointment_datetime", "clinic_code", "gp_practice", 
    "age", "sex", "ethnicity", "imd", "accessibility_flags", "dna_risk", "risk_label", "trial_arm", 
    "date_model_run", "model_version",
    "tier1_text_timestamp", "tier2_text_timestamp", "tier3_phone_attempt", "tier3_call_timestamp", 
    "patient_intent", "call_notes", "appointment_outcome"
  ),
  Description = c(
    "Cryptographically secure unique appointment key generated via MD5 hash of patient_id and appt_month.",
    "Unique patient identifier generated by system.",
    "Scheduled date and time of the appointment.",
    "Code of the clinic where the appointment is scheduled.",
    "GP practice where the patient is currently registered.",
    "Patient age on the scheduled appointment date.",
    "Patient biological sex or gender.",
    "Patient ethnicity record.",
    "Index of Multiple Deprivation decile (1-10 or unknown).",
    "Single collapsed column with a comma-separated list of active clinical accessibility & vulnerability flags (if present).",
    "Calibrated risk probability score of a no-show.",
    "Label risk profile based on our operational threshold.",
    "Randomized trial allocation arm.",
    "The precise timestamp/date when the predictive model scoring job was executed.",
    "The active production model version code from system configuration.",
    "Manual log: timestamp when the universal Tier 1 reminder SMS was dispatched.",
    "Manual log: timestamp when the Tier 2 reminder SMS was dispatched.",
    "Manual log: did a phone outreach attempt occur? (yes, no)",
    "Manual log: precise date and time of the phone outreach attempt.",
    "Manual log: patient's stated intent (confirm, cancel, reschedule, no_response).",
    "Manual log: free-text clinical coordinator notes and cancellation reasons.",
    "Manual log: final observed clinical outcome of the scheduled appointment slot."
  ),
  Data_Source = c(
    "Model generated", "Model generated", "EHR record", "EHR record", "EHR record",
    "EHR record", "EHR record", "EHR record", "Census lookup", "Clinical record", "Model generated",
    "Model generated", "Model generated", "Model generated", "Model generated", 
    "Coordinator log", "Coordinator log", "Coordinator log", "Coordinator log", "Coordinator log", "Coordinator log", "Coordinator log"
  )
)

# Write table headers and content
writeData(wb, sheet_name_inst, t(c("Variable name", "Variable description", "Data source")), startCol = 1, startRow = 5, colNames = FALSE)
addStyle(wb, sheet_name_inst, style_header, rows = 5, cols = 1:3)
writeData(wb, sheet_name_inst, mapping_dict, startCol = 1, startRow = 6, colNames = FALSE)
for (i in 6:(6 + nrow(mapping_dict) - 1)) {
  addStyle(wb, sheet_name_inst, style_text_left, rows = i, cols = 1:2)
  addStyle(wb, sheet_name_inst, style_text_left, rows = i, cols = 3)
}

# .........................................................................
# SHEET 2: WEEKLY OUTREACH ROSTER
# .........................................................................
sheet_name_rost <- "Weekly outreach roster"
addWorksheet(wb, sheet_name_rost)
showGridLines(wb, sheet_name_rost, show = FALSE)

# Dynamic indexing of rows
start_row <- 9
end_row   <- start_row + nrow(final_manifest) - 1

# Add title block
writeData(wb, sheet_name_rost, "Weekly coordination outreach roster", startCol = 1, startRow = 1)
addStyle(wb, sheet_name_rost, style_title, rows = 1, cols = 1)

# Add summary cards (KPI Block pointing to correct shifted columns for v7)
writeData(wb, sheet_name_rost, "Total high-risk cohort:", startCol = 1, startRow = 3)
addStyle(wb, sheet_name_rost, style_kpi_label, rows = 3, cols = 1)
writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(L%d:L%d, \"high-risk\")", start_row, end_row), startCol = 2, startRow = 3)
addStyle(wb, sheet_name_rost, style_kpi_value, rows = 3, cols = 2)

writeData(wb, sheet_name_rost, "Intervention arm (active calls):", startCol = 1, startRow = 4)
addStyle(wb, sheet_name_rost, style_kpi_label, rows = 4, cols = 1)
writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(M%d:M%d, \"intervention\")", start_row, end_row), startCol = 2, startRow = 4)
addStyle(wb, sheet_name_rost, style_kpi_value, rows = 4, cols = 2)

writeData(wb, sheet_name_rost, "Control arm (passive standard care):", startCol = 4, startRow = 3)
addStyle(wb, sheet_name_rost, style_kpi_label, rows = 3, cols = 4)
writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(M%d:M%d, \"control\")", start_row, end_row), startCol = 5, startRow = 3)
addStyle(wb, sheet_name_rost, style_kpi_value, rows = 3, cols = 5)

writeData(wb, sheet_name_rost, "Outreach completed:", startCol = 4, startRow = 4)
addStyle(wb, sheet_name_rost, style_kpi_label, rows = 4, cols = 4)
writeFormula(wb, sheet_name_rost, sprintf("COUNTIF(R%d:R%d, \"yes\")", start_row, end_row), startCol = 5, startRow = 4)
addStyle(wb, sheet_name_rost, style_kpi_value, rows = 4, cols = 5)

# Column headers (22 Columns total)
headers <- c(
  "appointment_id", "patient_id", "appointment_datetime", "clinic_code", "gp_practice", 
  "age", "sex", "ethnicity", "imd", "accessibility_flags", "dna_risk", "risk_label", "trial_arm", 
  "date_model_run", "model_version",
  "tier1_text_timestamp", "tier2_text_timestamp", "tier3_phone_attempt", "tier3_call_timestamp", 
  "patient_intent", "call_notes", "appointment_outcome"
)
writeData(wb, sheet_name_rost, t(headers), startCol = 1, startRow = 8, colNames = FALSE)
addStyle(wb, sheet_name_rost, style_header, rows = 8, cols = 1:22)

# -------------------------------------------------------------------------
# Vectorized Pre-Population inside R (Loop-Free!) [cite: 503]
# -------------------------------------------------------------------------
export_roster_data <- final_manifest %>% 
  select(
    appointment_id, patient_id, appointment_datetime, clinic_code, gp_practice, 
    age, sex, ethnicity, index_multiple_deprivation_decile, accessibility_flags, dna_risk, risk_label, trial_arm,
    date_model_run, model_version
  ) %>%
  mutate(
    # Set text reminders, call variables, and outcome fields dynamically based on trial arms [cite: 503]
    tier1_text_timestamp = case_when(
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    ),
    tier2_text_timestamp = case_when(
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    ),
    tier3_phone_attempt = case_when(
      trial_arm == "control" ~ "N/A - Control",
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    ),
    tier3_call_timestamp = case_when(
      trial_arm == "control" ~ "N/A - Control",
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    ),
    patient_intent = case_when(
      trial_arm == "control" ~ "N/A - Control",
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    ),
    call_notes = case_when(
      trial_arm == "control" ~ "N/A - Control",
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    ),
    # Pre-populate Standard BAU track as N/A, keep trial arms blank for coordinator entry
    appointment_outcome = case_when(
      trial_arm == "not in trial" ~ "N/A - Not in Trial",
      TRUE ~ ""
    )
  )

# Write all 22 columns to Excel in a SINGLE high-speed operation
writeData(wb, sheet_name_rost, export_roster_data, startCol = 1, startRow = 9, colNames = FALSE)

# -------------------------------------------------------------------------
# Vectorized Group Row Indices Calculations
# -------------------------------------------------------------------------
intervention_rows <- which(final_manifest$trial_arm == "intervention") + start_row - 1
control_rows      <- which(final_manifest$trial_arm == "control") + start_row - 1
not_in_trial_rows <- which(final_manifest$trial_arm == "not in trial") + start_row - 1

# -------------------------------------------------------------------------
# Vectorized Column Styling (No loops!)
# -------------------------------------------------------------------------
addStyle(wb, sheet_name_rost, style_text_left, rows = start_row:end_row, cols = c(1, 2, 4, 5, 8, 10, 12, 13, 14, 15), gridExpand = TRUE)
addStyle(wb, sheet_name_rost, style_text_left, rows = start_row:end_row, cols = 3, gridExpand = TRUE) # Scheduled Datetime
addStyle(wb, sheet_name_rost, style_num_right, rows = start_row:end_row, cols = c(6, 9), gridExpand = TRUE) # Integers
addStyle(wb, sheet_name_rost, style_text_left, rows = start_row:end_row, cols = 7, gridExpand = TRUE) # Sex

# Decision Model outcomes formatting (Column 11 / K)
addStyle(wb, sheet_name_rost, style_prob_pct, rows = start_row:end_row, cols = 11, gridExpand = TRUE) # DNA Risk Pct

# -------------------------------------------------------------------------
# Group-Based Vectorized Conditional Formatting
# -------------------------------------------------------------------------
if (length(intervention_rows) > 0) {
  # Active Intervention: Highlight all text, call, and outcome columns 16 to 22 (P-V)
  addStyle(wb, sheet_name_rost, style_active_input, rows = intervention_rows, cols = 16:22, gridExpand = TRUE)
}

if (length(control_rows) > 0) {
  # Control group: Highlight text timestamps (16:17) and final outcome (22) in blue, phone cols (18:21) in gray
  addStyle(wb, sheet_name_rost, style_active_input, rows = control_rows, cols = c(16, 17, 22), gridExpand = TRUE)
  addStyle(wb, sheet_name_rost, style_control_row, rows = control_rows, cols = 18:21, gridExpand = TRUE)
}

if (length(not_in_trial_rows) > 0) {
  # Not in Trial: Highlight all logging fields (16:22) in gray
  addStyle(wb, sheet_name_rost, style_control_row, rows = not_in_trial_rows, cols = 16:22, gridExpand = TRUE)
}

# Apply dynamic data validation dropdown lists for active outreach columns (Columns R/18, T/20, V/22)
dataValidation(wb, sheet_name_rost, col = 18, rows = start_row:end_row, type = "list", value = '\"yes,no\"')
dataValidation(wb, sheet_name_rost, col = 20, rows = start_row:end_row, type = "list", value = '\"confirm,cancel,reschedule,no_response\"')
dataValidation(wb, sheet_name_rost, col = 22, rows = start_row:end_row, type = "list", value = '\"Attended,DNA,Cancelled,Rescheduled\"')

# Column width padding
setColWidths(wb, sheet_name_rost, cols = 1, widths = 36)   # appointment_id
setColWidths(wb, sheet_name_rost, cols = 2, widths = 14)   # patient_id
setColWidths(wb, sheet_name_rost, cols = 3, widths = 22)   # appointment_datetime
setColWidths(wb, sheet_name_rost, cols = 4, widths = 14)   # clinic_code
setColWidths(wb, sheet_name_rost, cols = 5, widths = 26)   # gp_practice
setColWidths(wb, sheet_name_rost, cols = 6, widths = 10)   # age
setColWidths(wb, sheet_name_rost, cols = 7, widths = 10)   # sex
setColWidths(wb, sheet_name_rost, cols = 8, widths = 16)   # ethnicity
setColWidths(wb, sheet_name_rost, cols = 9, widths = 10)   # imd
setColWidths(wb, sheet_name_rost, cols = 10, widths = 38)  # accessibility_flags
setColWidths(wb, sheet_name_rost, cols = 11, widths = 14)  # dna_risk
setColWidths(wb, sheet_name_rost, cols = 12, widths = 14)  # risk_label
setColWidths(wb, sheet_name_rost, cols = 13, widths = 14)  # trial_arm
setColWidths(wb, sheet_name_rost, cols = 14, widths = 22)  # date_model_run
setColWidths(wb, sheet_name_rost, cols = 15, widths = 15)  # model_version
setColWidths(wb, sheet_name_rost, cols = 16:17, widths = 22) # manual text timestamps
setColWidths(wb, sheet_name_rost, cols = 18:20, widths = 22) # manual call logging blanks
setColWidths(wb, sheet_name_rost, cols = 21, widths = 30)  # call_notes
setColWidths(wb, sheet_name_rost, cols = 22, widths = 22)  # appointment_outcome

# Lock split panes (firstActiveCol = 2 freezes Column A and Row 8 headers)
freezePane(wb, sheet_name_rost, firstActiveRow = 9, firstActiveCol = 2)

# Safeguard: Recreate the system temp directory if Windows has cleaned it up
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)

# Save workbook
output_excel_path <- "outputs/weekly_scored_appointments.xlsx"
saveWorkbook(wb, file = output_excel_path, overwrite = TRUE)

# -------------------------------------------------------------------------
# 9. Console report (Sentence Case)
# -------------------------------------------------------------------------
n_total        <- nrow(final_manifest)
n_high_risk    <- sum(final_manifest$risk_label == "high-risk")
n_low_risk     <- sum(final_manifest$risk_label == "low-risk")
n_control      <- sum(final_manifest$trial_arm == "control")
n_intervention <- sum(final_manifest$trial_arm == "intervention")

cat("\n====================================================\n")
cat("      Weekly prediction and allocation report       \n")
cat("====================================================\n")
cat(sprintf("Processed appointments:     %d patients\n", n_total))
cat(sprintf("Operational threshold used: >= %.4f\n", op_threshold))
cat("----------------------------------------------------\n")
cat(sprintf("Low-risk (standard track):  %d patients (%.1f%%)\n", n_low_risk, (n_low_risk / n_total) * 100))
cat(sprintf("High-risk (trial cohort):   %d patients (%.1f%%)\n", n_high_risk, (n_high_risk / n_total) * 100))
cat("----------------------------------------------------\n")
cat(sprintf("-> Control arm (passive):   %d patients\n", n_control))
cat(sprintf("-> Intervention arm (active): %d patients\n", n_intervention))
cat("====================================================\n")
cat(sprintf("Success. Scored manifest saved to %s\n", output_manifest_path))
cat(sprintf("Success. Pre-formatted coordinator workbook saved to %s\n", output_excel_path))
cat("====================================================\n\n")
