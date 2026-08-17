source("00_libraries_and_utils.R")
source("02_data_prep.R")


if(tune_model) {
  source("03_tune_model.R")
} else {
  model_tune <- readRDS("data/processed/rf_tune.rds")
}

# fit best model according to PR

model <- local({
  
  params <- model_tune$best_params
  wf <- model_tune$tune_res %>% extract_workflow()
  fit <- wf %>%
    finalize_workflow(params) %>%
    fit(data = train_data)
  
  fit
})

vip::vip(model, 40)
