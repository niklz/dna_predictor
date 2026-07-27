
if(tune_model) {
  source("03_tune_model.R")
}

if(!tune_model) {
  source("02_data_prep.R")
  model_tune <- readRDS("data/rf_tune.RDS")
}

# fit best model according to PR

model <- local({
  
  params <- model_tune$best_params
  wf <- model_tune$fits %>% extract_workflow()
  fit <- wf %>%
    finalize_workflow(params) %>%
    fit(data = train_raw)
  
  fit
})

