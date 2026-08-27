# =========================================================================
# _MASTER_PREDICT.R: OPERATIONAL WEEKLY INFERENCE & ROSTER GENERATION 
# =========================================================================
# Orchestrates predictive scoring and workbook generation for both parent 
# organizations (SHOU and UHBW) using calibrated Random Forest models 

# 1. Setup and load configurations 
# -------------------------------------------------------------------------
source("00_libraries_and_utils.R")
conf <- config::get()

# Set seed for reproducible trial arm allocation 
set.seed(42) 

# 2. Load production artifacts 
# -------------------------------------------------------------------------
model_bundle_path <- "data/processed/rf_final_model.rds"
if (!file.exists(model_bundle_path)) {
  stop("Model bundle not found. Please train and calibrate the model first.")
}
model_bundle  <- readRDS(model_bundle_path)
rf_model      <- model_bundle$model
rf_calibrator <- model_bundle$calibrator

# Load calculated operational thresholds
threshold_path <- "data/processed/risk_threshold.RDS"
if (file.exists(threshold_path)) {
  threshold_data <- readRDS(threshold_path)
  # Standard default threshold
  calculated_threshold <- threshold_data$bounds$floor_threshold
} else {
  calculated_threshold <- 0.1446 # Power floor fallback if RDS is missing 
}

# Operational override: Manually setting threshold to 10% (0.1) to boost sample sizes
op_threshold <- 0.1 
message(sprintf("Target ML threshold locked at >= %.2f (Override active: %.2f calculated)", op_threshold, calculated_threshold))

# 3. Process SHOU Appointments [cite: 120]
# -------------------------------------------------------------------------
message("\n--- PROCESSING COHORT: SHOU ---")
raw_shou_path <- conf$new_appointment_path

if (!file.exists(raw_shou_path)) {
  # Fallback for dev/staging environments
  if (file.exists("data/data_joined.RDS")) {
    message("SHOU raw path not found. Initialising development fallback sample...")
    new_appts_shou <- read.csv("data/DNA_20260818.csv") %>% 
      sample_n(400) %>% 
      mutate(tumoursite = "Breast")
  } else {
    stop("SHOU raw appointments file and fallback database are both missing.")
  }
} else {
  new_appts_shou <- read.csv(raw_shou_path)
}

# Generate scored and randomized trial manifest 
final_manifest_shou <- generate_appointment_manifest(
  new_appointments = new_appts_shou,
  op_threshold     = op_threshold,
  rf_model         = rf_model,
  rf_calibrator    = rf_calibrator,
  conf             = conf
)

# Export formatted Excel workbook and scored CSV lists for SHOU 
generate_excel_manifest(
  manifest     = final_manifest_shou, 
  org_name     = "shou", 
  op_threshold = op_threshold
)


# 4. Process UHBW Appointments 
# -------------------------------------------------------------------------
message("\n--- PROCESSING COHORT: UHBW ---")
raw_uhbw_path <- "data/new_clinic_appointments_uhbw.xlsx"

if (!file.exists(raw_uhbw_path)) {
  stop("UHBW Raw appointment file not found at: ", raw_uhbw_path)
}

# Cleanly read Excel file and force all incoming headers to lowercase
new_appts_uhbw <- readxl::read_excel(raw_uhbw_path) %>% 
  rename_with(.fn = str_to_lower)

# Generate scored and randomized trial manifest 
final_manifest_uhbw <- generate_appointment_manifest(
  new_appointments = new_appts_uhbw,
  op_threshold     = op_threshold,
  rf_model         = rf_model,
  rf_calibrator    = rf_calibrator,
  conf             = conf
)

# Export formatted Excel workbook and scored CSV lists for UHBW 
generate_excel_manifest(
  manifest     = final_manifest_uhbw, 
  org_name     = "uhbw", 
  op_threshold = op_threshold
)

message("\nInference Pipeline Completed Successfully! Scored weekly outputs are ready for coordinator distribution.")
