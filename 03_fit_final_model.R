source("00_libraries_and_utils.R")
source("01_data_prep.R")
conf <- config::get()

train_data <- readRDS("data/processed/train_engineered.rds")


if(conf$tune_model) {
  source("02_tune_model.R")
} else {
  model_tune_results <- readRDS("data/processed/rf_tune.rds")
}

# fit best model according to PR

model <- local({
  
  params <- model_tune_results$best_params
  wf <- model_tune_results$tune_res %>% extract_workflow()
  
  # ---------------------------------------------------------
  # 1. Speed up Final Model Training (Multi-threading Override)
  # ---------------------------------------------------------
  # Extract the current parsnip spec and update ranger's threads [cite: 419]
  current_spec <- extract_spec_parsnip(wf)
  
  updated_spec <- current_spec %>%
    set_engine(
      "ranger",
      # Dynamically detects cores and uses all but one for OS stability
      num.threads = parallel::detectCores() - 1, 
      importance = "permutation"
    )
  
  wf <- wf %>% 
    remove_model() %>% 
    add_model(updated_spec)
  
  fit <- wf %>%
    finalize_workflow(params) %>%
    fit(data = train_data)

  preds <- model_tune_results$tune_res %>%
    collect_predictions(parameters  = params)
  
  cal_model <- cal_estimate_logistic(
    preds,
    truth = dna_outcome,
    estimate = c(.pred_DNA, .pred_attended),
    event_level = "first"
  )
  
  
  list(
    model = fit,
    calibrator = cal_model,
    model_ver = conf$model_ver
  )
})

saveRDS(model, "data/processed/rf_final_model.rds")
