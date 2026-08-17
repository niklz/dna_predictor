dataset <- read.csv("data/dna_data.csv")

ds <- dataset[dataset$Historic == 1, ]

ds$dna_outcome <- ifelse(ds$AttendedStatusCode == "3", 1, 0)
ds$prev_dna <- ifelse(is.na(ds$prev_dna), 0, ds$prev_dna)
ds$imd <- as.character(ds$Index_Multiple_Deprivation_Decile)
ds$imd <- ifelse(is.na(ds$imd), "unknown", ds$imd)

# split data
set.seed(123)

rows <- nrow(ds)
prop <- 0.7
idx <- sample(seq_len(nrow(ds)), size = floor(rows*prop))
ds_train <- ds[idx, ]
ds_test <- ds[-idx, ]

model <- glm(dna_outcome ~ AgeAtAppointment + Gender + Autistic.spectrum.disorder + imd + prev_dna + ConsultationMedia, data = ds_train, family = binomial)

pred <- predict(model, newdata = ds_test,  type = "response")

actual <- ds_test$dna_outcome

roc_data <- data.frame(actual, pred)
roc_data <- roc_data[order(-roc_data$pred), ]

# Calculate TPR and FPR at each threshold
thresholds <- sort(unique(roc_data$pred), decreasing = TRUE)
tpr <- numeric(length(thresholds))
fpr <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
  cutoff <- thresholds[i]
  predicted_class <- ifelse(roc_data$pred >= cutoff, 1, 0)
  
  tp <- sum(predicted_class == 1 & roc_data$actual == 1)
  fp <- sum(predicted_class == 1 & roc_data$actual == 0)
  fn <- sum(predicted_class == 0 & roc_data$actual == 1)
  tn <- sum(predicted_class == 0 & roc_data$actual == 0)
  
  tpr[i] <- tp / (tp + fn)  # Sensitivity
  fpr[i] <- fp / (fp + tn)  # 1 - Specificity
}

# Plot ROC curve
plot(fpr, tpr, type = "l", col = "blue", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "ROC Curve")
abline(0, 1, col = "gray", lty = 2)  # Diagonal reference line
auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)


ds$prob_dna <- pred 

dataset <- ds
