source("00_libraries_and_utils.R")
source("01_data_prep.R")
conf <- config::get()


if(conf$tune_model) {
  source("02_tune_model.R")
} else {
  model_tune_results <- readRDS("data/processed/rf_tune.rds")
}

# fit best model according to PR

model <- local({
  
  params <- model_tune_results$best_params
  wf <- model_tune_results$tune_res %>% extract_workflow()
  fit <- wf %>%
    finalize_workflow(params) %>%
    fit(data = train_data)
  
  fit
})


saveRDS(model, "data/processed/rf_final_model.rds")
