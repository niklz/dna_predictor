source("00_libraries_and_utils.R")

# -------------------------------------------------------------------------
# Remote Progress Monitoring: Mirror console to the network drive
# -------------------------------------------------------------------------
network_log_dir <- "S:/Finance/Shared Area/BNSSG - BI/8 Modelling and Analytics/working/nh/projects/dna_predictor/data"
log_file_path   <- file.path(network_log_dir, "live_training_progress.log")

# Create the file connection
log_connection <- file(log_file_path, open = "wt")

# Sink console output and messages to the network log file [cite: 7]
sink(log_connection, split = TRUE)          # split = TRUE keeps output visible in RStudio!
sink(log_connection, type = "message")       # Redirects messages and warnings safely


conf <- config::get()
# source("01_data_prep.R")
source("02_tune_model.R")
source("03_fit_final_model.R")
source("04_set_threshold.R")

# Clean up and close console redirection
sink(type = "message")
sink()
close(log_connection)
message("Run complete. Console logging dismissed.")
