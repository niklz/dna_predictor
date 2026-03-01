source("02_data_prep.R")

mod_out <- local({

mod_test <- model_data$test
mod_train <- model_data$train

target <- "dna_outcome"
features <- names(mod_train)[! names(mod_train) %in% target]

# quick and dirty feature importance ranking

results <- data.frame(
  feature = features,
  form = NA,
  performance = NA,
  status = "Pending",
  stringsAsFactors = FALSE
) 

train_samp <- mod_train[sample(seq_len(nrow(mod_train)), size = n_featimp ), ] 

for (i in seq_along(features)) {
  current_feature <- features[i]
  
  form_str <- paste(target, "~", current_feature)
  form <- as.formula(form_str)
  
  results$form[i] <- form_str
  
  tryCatch({
    # Run the model
    mod <- glm(form, data = train_samp, family = binomial)
    
    results$performance[i] <- summary(mod)$aic
    results$status[i]      <- "Success"
    
  }, error = function(e) {
    results$status[i]      <- paste("Error:", conditionMessage(e))
  })
}


top_feat <- results[order(results$performance), ]$feature

form <- as.formula(paste(target, "~", paste(top_feat[-1], collapse = " + ")))

# train model
model <- glm(form, data = mod_train, family = binomial)

# generate predictions
pred <- predict(model, newdata = mod_test,  type = "response")
actual <- mod_test$dna_outcome

# ROC curve
roc_data <- data.frame(actual = mod_test$dna_outcome, pred = pred)
roc_data <- roc_data[order(-roc_data$pred), ]

tp_cum <- cumsum(roc_data$actual == 1)
fp_cum <- cumsum(roc_data$actual == 0)

total_p <- sum(roc_data$actual == 1)
total_n <- sum(roc_data$actual == 0)

tpr <- tp_cum / total_p
fpr <- fp_cum / total_n

roc_plot <- plot(fpr, tpr, type = "l", col = "blue", lwd = 2,
                 xlab = "False Positive Rate", ylab = "True Positive Rate",
                 main = "Optimized ROC Curve")
abline(0, 1, col = "gray", lty = 2)

auc <- sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2)

# importance
model_summary <- summary(model)
coeffs <- model_summary$coefficients
importance <- abs(coeffs[, 3])
importance <- importance[names(importance) != "(Intercept)"]
sorted_importance <- sort(importance, decreasing = TRUE)
print(sorted_importance)

feat_imp_plot <- barplot(sorted_importance, 
        las = 2, 
        col = "skyblue", 
        main = "Feature Importance (Absolute z-score)",
        ylab = "Absolute z-value")

list(model = model, auc = auc, roc_plot = roc_plot, feat_imp_plot = feat_imp_plot)

})
