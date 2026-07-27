## hyper params for modelling

# # number of levels to lump factors into
# n_lump <- 10
# # test/train split
# prop_tt <- 0.7
# # quick feature importance rows to train
# n_featimp <- 5E3
# # near zero variance frequency
# nzv_freq <- 0.80

# model specs / vars

target_col <- "dna_outcome"

vars <- c(
  "local_spec_code",
  "national_spec_code",
  "appointment_type",
  "distance_km",
  "nfa_ind",
  "age_group",
  "age_at_appointment",
  "ethnicity",
  "a_ld",
  "a_autism",
  "a_interpreter_req_bsl",
  "a_interpreter_req_lang",
  "a_balance",
  "a_cognitive_impairment",
  "a_mobility_restriction",
  "a_hear_vis_impaired",
  "a_dementia",
  "a_depression",
  "a_downs_syndrome",
  "a_long_standing_condition",
  "a_makaton",
  "a_mild_cognitive_impairment",
  "a_memory_impairment",
  "a_mood_disorder",
  "a_other_disability",
  "a_psychosis",
  "a_severe_anxiety",
  "a_wheelchair_user",
  "gender",
  "registered_gp_practice",
  "site_code",
  "prev_dna_ly",
  "appt_hour",
  "appt_dow",
  "appt_month",
  "appt_wknd_ind",
  "referral_urgency",
  "lead_time_days",
  "clinic_code",
  "clinic_location",
  "imd"
)

# Flag to tune model params over grid (this is slow), otherwise use pre-trained model
tune_model <- FALSE
# Number of worker to parallelise (note that ranger already uses 4 threads per worker so set this accordingly!)
num_workers <- 6

# RF hyperparams
tree_range <- c(250, 1000)
mtry <- 2
trees <- 468
min_n <- 6
fct_other_prp <- 0.02

model_ver <- "DEV"