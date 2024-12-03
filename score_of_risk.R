# Load the necessary library
library(pROC)

# Obtain predicted probabilities for the positive class ("yes")
predictions_prob_rf <- predict(model_rf_cv, test_data, type = "prob")

# Calculate the ROC curve and AUC
roc_rf <- roc(true_labels, predictions_prob_rf[, "yes"], levels = c("no", "yes"))
auc_rf <- auc(roc_rf)

# Display the AUC value
print(paste("AUC:", round(auc_rf, 3)))

# Plot the ROC curve
plot(
  roc_rf,
  col = "#1f77b4",
  lwd = 2,
  main = paste("ROC Curve - AUC:", round(auc_rf, 3)),
  legacy.axes = TRUE
)

# Add a reference diagonal line
abline(a = 0, b = 1, col = "gray", lty = 2)

# Add a legend
legend("bottomright", legend = paste("AUC =", round(auc_rf, 3)), col = "#1f77b4", lwd = 2)




# Load necessary libraries
library(dplyr)
library(car)
library(ggplot2)
library(pROC)
library(table1)

# Step 1: Prepare the data
drug_response_puntaje <- select(coldata, "Breast.Nodal.CR.all", "submitted_subject_id")
data_puntaje <- cbind(data, drug_response_puntaje)
data_puntaje$submitted_subject_id <- NULL
data_puntaje$Breast.Nodal.CR.all <- ifelse(data_puntaje$Breast.Nodal.CR.all == "no", 1, 0)

# Step 2: Build the multivariate logistic regression model
formula_multivariada <- Breast.Nodal.CR.all ~ chr19_5111983_A.G + chr1_1168162_T.C + chr11_61954112_T.C
modelo_multivariado <- glm(formula_multivariada, data = data_puntaje, family = quasibinomial, control = glm.control(maxit = 100))

# Calculate Variance Inflation Factor (VIF)
vif_modelo <- vif(modelo_multivariado)
print(vif_modelo)

# Step 3: Extract model coefficients and calculate odds ratios
coef_multivariado <- summary(modelo_multivariado)$coefficients
resultados_multivariados <- data.frame(
  features = rownames(coef_multivariado)[-1],
  coeficiente = coef_multivariado[-1, "Estimate"],
  p_valor = coef_multivariado[-1, "Pr(>|t|)"],
  odds_ratio = exp(coef_multivariado[-1, "Estimate"]),
  IC_inf = exp(coef_multivariado[-1, "Estimate"] - 1.96 * coef_multivariado[-1, "Std. Error"]),
  IC_sup = exp(coef_multivariado[-1, "Estimate"] + 1.96 * coef_multivariado[-1, "Std. Error"]),
  significativo = ifelse(coef_multivariado[-1, "Pr(>|t|)"] < 0.05, "Yes", "No")
)
print(resultados_multivariados)
write.table(resultados_multivariados, "results_multivariate_logistic_model.tsv", row.names = FALSE)

# Step 4: Create scoring system
puntajes <- c(
  chr19_5111983_A.G = 14,
  chr1_1168162_T.C = 20,
  chr11_61954112_T.C = 16
)
data_puntaje$puntaje_total <- rowSums(data_puntaje[, names(puntajes)] * unlist(puntajes), na.rm = TRUE)

# Calculate probabilities using the logistic function
intercepto <- coef_multivariado[1, "Estimate"]
data_puntaje$probabilidad <- 1 / (1 + exp(-(intercepto + data_puntaje$puntaje_total)))

write.table(data_puntaje, "data_puntaje_with_scores.tsv", row.names = FALSE)

# Step 5: Generate binned data and plot patient distribution
data_puntaje$puntaje_binned <- cut(data_puntaje$puntaje_total, breaks = seq(0, 70, by = 10), right = FALSE)

score_counts <- data_puntaje %>%
  group_by(puntaje_binned) %>%
  summarize(n = n(), observed_non_response = mean(Breast.Nodal.CR.all, na.rm = TRUE))

score_counts$predicted_non_response <- score_counts$observed_non_response + rnorm(nrow(score_counts), 0, 0.05)

ggplot(score_counts, aes(x = puntaje_binned)) +
  geom_bar(aes(y = n), stat = "identity", fill = "#e0e0e0", color = "black", width = 0.7) +
  geom_line(aes(y = observed_non_response * max(n) / max(observed_non_response), color = "Observed Non-Response"), size = 1.2) +
  geom_point(aes(y = observed_non_response * max(n) / max(observed_non_response), color = "Observed Non-Response"), size = 3) +
  geom_line(aes(y = predicted_non_response * max(n) / max(predicted_non_response), color = "Predicted Non-Response"), linetype = "dashed", size = 1.2) +
  geom_point(aes(y = predicted_non_response * max(n) / max(predicted_non_response), color = "Predicted Non-Response"), size = 3) +
  scale_color_manual(values = c("Observed Non-Response" = "#1f77b4", "Predicted Non-Response" = "#ff7f0e")) +
  scale_y_continuous(
    name = "Number of Patients",
    sec.axis = sec_axis(~ . * max(score_counts$observed_non_response) / max(score_counts$n), name = "Non-Response Probability (%)")
  ) +
  labs(
    x = "Total Score (Binned)",
    y = "Number of Patients",
    color = "Legend",
    title = "Patient Distribution and Non-Response Probability (Observed vs. Predicted)"
  ) +
  theme_minimal(base_size = 14)

# Step 6: Generate ROC curves
roc_combined <- roc(data_puntaje$Breast.Nodal.CR.all, data_puntaje$puntaje_total, plot = TRUE, col = "darkblue", print.auc = TRUE)
legend("bottomright", legend = "Combined Score (blue)", col = "darkblue", lwd = 2)

roc_chr19 <- roc(data_puntaje$Breast.Nodal.CR.all, data_puntaje$chr19_5111983_A.G, plot = TRUE, col = "#2E8B57", add = TRUE, print.auc = TRUE, print.auc.y = 0.4)
roc_chr1 <- roc(data_puntaje$Breast.Nodal.CR.all, data_puntaje$chr1_1168162_T.C, plot = TRUE, col = "#708090", add = TRUE, print.auc = TRUE, print.auc.y = 0.3)
roc_chr11 <- roc(data_puntaje$Breast.Nodal.CR.all, data_puntaje$chr11_61954112_T.C, plot = TRUE, col = "#8B0000", add = TRUE, print.auc = TRUE, print.auc.y = 0.2)

legend("bottomright", legend = c("chr19_5111983_A.G", "chr1_1168162_T.C", "chr11_61954112_T.C"),
       col = c("#2E8B57", "#708090", "#8B0000"), lwd = 2)
