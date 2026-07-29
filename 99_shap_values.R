library(fastshap)
library(shapviz)
library(tidymodels)
library(doParallel)
library(parallel)

# 1. Setup Parallel Cluster
# Detect cores and leave 1 core free for system background tasks
n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores)

# Ensure worker processes can find installed libraries (important on Windows)
invisible(clusterCall(cl, function(lp) .libPaths(lp), .libPaths()))
registerDoParallel(cl)

# 1. Extract raw ranger engine & recipe from your existing model
raw_rf <- extract_fit_engine(model)
rec    <- extract_recipe(model)

# 2. Bake your data (features only, target removed)
train_baked <- bake(rec, new_data = train_raw) %>% 
  select(-any_of(target_col)) # Replace target_col with your target column name

# 3. Background dataset (Sample 300–500 rows for fast background sampling)
bg_X <- train_baked[sample(nrow(train_baked), 300), ]

pfun <- function(object, newdata) {
  # Extract predictions for the positive target class (Column 2)
  predict(object, data = newdata)$predictions[, 2]
}

# Pick 1 observation to explain (e.g., the first row of your dataset)
single_row <- train_baked[7483, , drop = FALSE]

# Fast local explanation (Completes in < 1 second!)
ex_local <- explain(
  object = raw_rf,
  X = bg_X,                      # Background comparison data
  newdata = single_row,          # THE KEY ARGUMENT: Only explain this single row
  pred_wrapper = pfun,
  nsim = 100,
  adjust = TRUE                  # Fixes sum efficiency
)

# Plot Waterfall using shapviz
baseline <- mean(pfun(raw_rf, newdata = bg_X))
shv <- shapviz(ex_local, X = single_row, baseline = baseline)

sv_waterfall(shv)

# Sample 1000 rows to explain for the global plot
explain_sample <- train_baked[sample(nrow(train_baked), 1000), ]

# 3. Run fastshap in parallel
# Note the two key parameters: parallel = TRUE and .packages = "ranger"
ex_global <- explain(
  object = raw_rf,
  X = bg_X,                      # Background matrix (e.g., 300 rows)
  newdata = explain_sample,      # Target rows to explain (e.g., 1000 rows)
  pred_wrapper = pfun,
  nsim = 50,                     # Higher simulation count for smoother results
  adjust = TRUE,
  parallel = TRUE,               # <--- ENABLES PARALLELISM
  .packages = "ranger"           # <--- CRITICAL: Loads ranger onto worker nodes
)

# 4. Stop Cluster (Always clean up background processes!)
stopCluster(cl)
registerDoSEQ()                  # Reset R back to standard single-threaded mode

# 5. Visualize
shv_global <- shapviz(ex_global, X = explain_sample)
sv_importance(shv_global, kind = "beeswarm")

# Convert and plot Beeswarm
shv_global <- shapviz(ex_global, X = explain_sample)
sv_importance(shv_global, kind = "beeswarm")
