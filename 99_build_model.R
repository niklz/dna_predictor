source("02_data_prep.R")

mod_out <- local({
  mod_test <- model_data$test
  mod_train <- model_data$train

  mod_train$dna_outcome <- factor(mod_train$dna_outcome, levels = c("1", "0"))
  mod_test$dna_outcome <- factor(mod_test$dna_outcome, levels = c("1", "0"))

  target <- "dna_outcome"
  features <- names(mod_train)[!names(mod_train) %in% target]

  # quick and dirty feature importance ranking

  results <- data.frame(
    feature = features,
    form = NA,
    performance = NA,
    status = "Pending",
    stringsAsFactors = FALSE
  )

  train_samp <- mod_train[sample(seq_len(nrow(mod_train)), size = n_featimp), ]

  for (i in seq_along(features)) {
    current_feature <- features[i]

    form_str <- paste(target, "~", current_feature)
    form <- as.formula(form_str)

    results$form[i] <- form_str

    tryCatch(
      {
        # Run the model
        mod <- glm(form, data = train_samp, family = binomial)

        results$performance[i] <- summary(mod)$aic
        results$status[i] <- "Success"
      },
      error = function(e) {
        results$status[i] <- paste("Error:", conditionMessage(e))
      }
    )
  }

  top_feat <- results[order(results$performance), ]$feature

  form <- as.formula(paste(target, "~", paste(top_feat[-1], collapse = " + ")))

  model <- ranger(
    form,
    data = mod_train,
    importance = 'permutation',
    sample.fraction = c(0.5, 0.5),
    replace = TRUE,
    probability = TRUE,
    mtry = mtry, # Try more features per split
    min.node.size = nodes # Prevents overfitting on noise
  )

  # generate predictions
  pred <- predict(model, data = mod_test, type = "response")
  prob_dna <- pred$predictions[, "1"]
  actual <- as.numeric(as.character(mod_test$dna_outcome))

  # ROC curve
  roc_data <- data.frame(actual = actual, pred = prob_dna)
  roc_data <- roc_data[order(-roc_data$pred), ]

  tp_cum <- cumsum(roc_data$actual == 1)
  fp_cum <- cumsum(roc_data$actual == 0)

  total_p <- sum(roc_data$actual == 1)
  total_n <- sum(roc_data$actual == 0)

  tpr <- tp_cum / total_p
  fpr <- fp_cum / total_n

  roc_plot <- ggplot(data.frame(fpr = fpr, tpr = tpr), aes(x = fpr, y = tpr)) +
    geom_line(color = "blue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = "Optimised ROC Curve",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal() +
    coord_fixed()
  roc_plot

  auc <- sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2)

  # importance
  importance <- vi(model)
  feat_imp_plot <- vip(importance)

  list(
    model = model,
    auc = auc,
    roc_plot = roc_plot,
    feat_imp_plot = feat_imp_plot
  )
})
