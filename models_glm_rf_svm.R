###### Required Packages ########
library(yardstick)
library(caret)
library(MLmetrics)
library(boot)
library(dplyr)
library(glmnet)
library(randomForest)
library(FactoMineR)
library(ggplot2)

######### Data Preparation and integration #########

# Merge dataframe of each omics: deg and/or tDNA and/or RNA editing -> "merged_data"


##### Clinical data preparation ####################

# Convert to dummy variables (one-hot encoding)
dummy_data <- dummyVars("~ . - submitted_subject_id - Breast.Nodal.CR.all", data = selection)
data_clinical <- as.data.frame(predict(dummy_data, newdata = selection))
data_clinical$submitted_subject_id <- selection$submitted_subject_id



######## PCA ########

gene_ids <- merged_data$gene_id
merged_data_no_gene_id <- merged_data[, -which(colnames(merged_data) == "gene_id")]
transposed_countdata <- t(merged_data_no_gene_id)
colnames(transposed_countdata) <- gene_ids
variances <- apply(transposed_countdata, 2, var)
filtered_transposed_countdata <- transposed_countdata[, variances != 0]
pca_result <- prcomp(filtered_transposed_countdata, scale = TRUE)

# Plot PCA results
plot(cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2)), type = "b", 
     xlab = "Principal Components", ylab = "Cumulative Variance Explained")

######## LASSO Regression ########

response_variable <- selection$Breast.Nodal.CR.all
response_variable_numeric <- ifelse(response_variable == "no", 1, 0)
merged_transposed_df <- as.data.frame(t(merged_data_no_gene_id))
lasso_model <- cv.glmnet(as.matrix(merged_transposed_df), response_variable_numeric, alpha = 1)
important_variables <- as.data.frame(as.matrix(coef(lasso_model, s = "lambda.min"))[coef(lasso_model, s = "lambda.min") != 0, ])
filtered_data <- merged_transposed_df[, rownames(important_variables)]

######## GLM Regression ########

# Create logistic regression model
train_indices <- sample(nrow(filtered_data), 0.7 * nrow(filtered_data))
train_data <- filtered_data[train_indices, ]
test_data <- filtered_data[-train_indices, ]
train_labels <- response_variable_numeric[train_indices]
test_labels <- response_variable_numeric[-train_indices]

train_labels <- as.factor(train_labels)

model_glm_cv <- train(train_labels ~ ., 
                      data = train_data, 
                      method = "glmnet", 
                      family = "binomial", 
                      trControl = trainControl(method = "cv", number = 10),
                      tuneGrid = expand.grid(alpha = 0:1, lambda = seq(0.001, 0.1, length = 10)))

# Make predictions on test set
predictions_glm <- predict(model_glm_cv, test_data)

# Evaluate model performance
df_glm <- data.frame(Prediction = as.numeric(predictions_glm), True_Labels = as.numeric(test_labels))
df_glm$Prediction <- factor(df_glm$Prediction, levels = c(0, 1))
df_glm$True_Labels <- factor(df_glm$True_Labels, levels = c(0, 1))

# Compute F1 score
f1_glm <- F1_Score(df_glm$Prediction, df_glm$True_Labels)

######## Random Forest and SVM ########

# Train Random Forest
model_rf_cv <- train(response_variable ~ ., data = filtered_data, method = "rf", trControl = trainControl(method = "cv", number = 10), importance = TRUE)

# Train SVM
model_svm_cv <- train(response_variable ~ ., data = filtered_data, method = "svmLinear", trControl = trainControl(method = "cv", number = 10))

######## Bootstrap for F1-Score Confidence Intervals ########

# Function to calculate F1 score in bootstrap samples
bootstrap_f1 <- function(data, indices) {
  sample_data <- data[indices, ]
  F1_Score(sample_data$Prediction, sample_data$True_Labels)
}

# Logistic Regression
bootstrap_results_glm <- boot(data = df_glm, statistic = bootstrap_f1, R = 1000)
f1_scores_glm <- bootstrap_results_glm$t
ci_glm <- boot.ci(bootstrap_results_glm, type = "basic")

# Random Forest
df_rf <- data.frame(Prediction = predict(model_rf_cv, test_data), True_Labels = test_labels)
bootstrap_results_rf <- boot(data = df_rf, statistic = bootstrap_f1, R = 1000)
f1_scores_rf <- bootstrap_results_rf$t
ci_rf <- boot.ci(bootstrap_results_rf, type = "basic")

# SVM
df_svm <- data.frame(Prediction = predict(model_svm_cv, test_data), True_Labels = test_labels)
bootstrap_results_svm <- boot(data = df_svm, statistic = bootstrap_f1, R = 1000)
f1_scores_svm <- bootstrap_results_svm$t
ci_svm <- boot.ci(bootstrap_results_svm, type = "basic")
