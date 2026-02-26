##############################################################################
# Multi-omics ML pipeline (clean version for GitHub)
# - Same steps as your original script (no methodological changes)
# - Organized, deduplicated, and made paths GitHub-friendly
# - Tip: put input files under ./data/... and results under ./results/...
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(glmnet)
  library(randomForest)
  library(kernlab)     # svmLinear
  library(MLmetrics)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(progress)
  library(patchwork)
  library(grid)
})

# Optional (used later). Install if missing.
if (!requireNamespace("PRROC", quietly = TRUE)) {
  install.packages("PRROC")
}

##############################################################################
# 0) PATHS (replace with your repo structure)
##############################################################################
# Prefer relative paths for GitHub (avoid ~/Desktop/... and setwd()).
# Example layout:
#   data/DEG/all_data_pivoted_104_DEG.tsv
#   data/clinical/clinical_data_forR.tsv
#   data/clinical/BEAUTY_SraRunTable.txt
#   data/clinical/phs001050.v1.p1.txt
#   data/data_analysis/tdna_gdna_290_47_rn.tsv
#   data/reanalisis/reditool_final.csv
#   results/

path_deg       <- "data/DEG/all_data_pivoted_104_DEG.tsv"
path_clin      <- "data/clinical/clinical_data_forR.tsv"
path_sra       <- "data/clinical/BEAUTY_SraRunTable.txt"
path_phs       <- "data/clinical/phs001050.v1.p1.txt"      # loaded (kept, as in your script)
path_dna       <- "data/data_analysis/tdna_gdna_290_47_rn.tsv"
path_reditool  <- "data/reanalisis/reditool_final.csv"

dir.create("results", showWarnings = FALSE, recursive = TRUE)

##############################################################################
# 1) LOAD INPUTS
##############################################################################
deg_raw <- read.csv(path_deg, sep = "")
rownames(deg_raw) <- deg_raw$Name
deg_raw$Name <- NULL
Run1 <- colnames(deg_raw)

clinical_data_forR <- read.delim(path_clin)
BEAUTY_SraRunTable  <- read.csv(path_sra)
phs001050.v1.p1     <- read.delim(path_phs)   # kept (even if not used later)

##############################################################################
# 2) MERGE CLINICAL + FILTER SAMPLES
##############################################################################
merge_tbl <- merge(
  BEAUTY_SraRunTable, clinical_data_forR,
  by.x = c("submitted_subject_id"),
  by.y = c("Subject_ID"),
  all = TRUE
)

merge_tbl <- merge_tbl %>% filter(Run %in% Run1)
merge_filtered <- merge_tbl %>% filter(analyte_type == "RNA")

selection <- dplyr::select(
  merge_filtered,
  submitted_subject_id, histological_type, Breast.Nodal.CR.all,
  Clinical.Molecular.Subtype, T.stage, N.stage, Age.group
)

# Dummy variables (one-hot), excluding ID and outcome
dummy_data <- dummyVars("~ . - submitted_subject_id - Breast.Nodal.CR.all", data = selection)
data_transformed <- predict(dummy_data, newdata = selection)
data_clinical <- as.data.frame(data_transformed)
data_clinical$submitted_subject_id <- selection$submitted_subject_id

##############################################################################
# 3) LOAD OMICS + STANDARDIZE ID COLUMN
##############################################################################
deg <- read.csv(path_deg, sep = "")
dna <- read.csv(path_dna, sep = "")
reditool <- read.csv(path_reditool, sep = "")

deg <- deg %>% dplyr::rename(gene_id = Name)
reditool <- reditool %>% dplyr::rename(gene_id = Uploaded_variation)

# (Optional checks retained)
length(unique(data_clinical$submitted_subject_id))
colnames(deg)[1:5]; head(deg[, 1:3])
colnames(reditool)[1:5]; head(reditool[, 1:3])
colnames(dna)[1:5]; head(dna[, 1:3])
"gene_id" %in% colnames(deg)
"gene_id" %in% colnames(reditool)
"gene_id" %in% colnames(dna)

##############################################################################
# 4) MAP Run (SRR) -> submitted_subject_id AND RENAME OMICS COLUMNS
##############################################################################
map <- merge_filtered %>%
  dplyr::select(submitted_subject_id, Run) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(Run), !is.na(submitted_subject_id))

srr_to_ex <- setNames(map$submitted_subject_id, map$Run)

rename_cols_by_map <- function(df, srr_to_ex) {
  old <- colnames(df)
  srr_cols <- intersect(old, names(srr_to_ex))
  if (length(srr_cols) > 0) {
    colnames(df)[match(srr_cols, colnames(df))] <- unname(srr_to_ex[srr_cols])
  }
  df
}

deg <- rename_cols_by_map(deg, srr_to_ex)
reditool <- rename_cols_by_map(reditool, srr_to_ex)

# Align IDs across clinical + 3 omics
ids_ok <- Reduce(intersect, list(
  data_clinical$submitted_subject_id,
  colnames(deg)[-1],
  colnames(reditool)[-1],
  colnames(dna)[-1]
))

deg <- deg[, c("gene_id", ids_ok), drop = FALSE]
reditool <- reditool[, c("gene_id", ids_ok), drop = FALSE]
dna <- dna[, c("gene_id", ids_ok), drop = FALSE]

data_clinical <- data_clinical %>% dplyr::filter(submitted_subject_id %in% ids_ok)
selection <- selection %>% dplyr::filter(submitted_subject_id %in% ids_ok)

# Prefixes (same as your script)
deg$gene_id <- paste0("EXP__", deg$gene_id)
reditool$gene_id <- paste0("ED__", reditool$gene_id)
# dna keeps its existing gene_id prefixes (gDNA_/tDNA_), as in your notes

merged_list <- list(
  EXP = deg,
  ED = reditool,
  DNA = dna,
  ED_EXP = dplyr::bind_rows(reditool, deg),
  ED_DNA = dplyr::bind_rows(reditool, dna),
  EXP_DNA = dplyr::bind_rows(deg, dna),
  ED_EXP_DNA = dplyr::bind_rows(reditool, deg, dna)
)

head(names(srr_to_ex))
head(colnames(deg)[1:6])
head(colnames(reditool)[1:6])

##############################################################################
# 5) HELPERS (same logic)
##############################################################################
safe_p <- function(x, y_bin01) {
  if (length(unique(y_bin01)) < 2) return(1)
  if (sd(x, na.rm = TRUE) == 0) return(1)
  tryCatch(t.test(x ~ y_bin01)$p.value, error = function(e) 1)
}

build_X_omics <- function(merged_data) {
  stopifnot("gene_id" %in% colnames(merged_data))
  gene_ids <- merged_data$gene_id
  X <- t(merged_data[, -1, drop = FALSE])
  colnames(X) <- gene_ids
  rownames(X) <- colnames(merged_data)[-1]
  X <- as.data.frame(X)
  rn <- rownames(X)
  X[] <- lapply(X, function(z) as.numeric(as.character(z)))
  rownames(X) <- rn
  X[is.na(X)] <- 0
  X
}

align_clinical <- function(X, data_clinical) {
  stopifnot("submitted_subject_id" %in% colnames(data_clinical))
  dc <- data_clinical
  dc$submitted_subject_id <- as.character(dc$submitted_subject_id)
  dc2 <- dc[match(rownames(X), dc$submitted_subject_id), , drop = FALSE]
  stopifnot(identical(rownames(X), dc2$submitted_subject_id))
  dc2
}

align_y <- function(X, data_clinical, y_vec) {
  y <- factor(y_vec, levels = c("yes", "no"))
  names(y) <- as.character(data_clinical$submitted_subject_id)
  y2 <- y[rownames(X)]
  stopifnot(!any(is.na(y2)))
  y2
}

##############################################################################
# 6) PIPELINE (same steps: split -> var0 train -> diff train -> lasso train -> CV models -> test)
##############################################################################
run_pipeline <- function(merged_data, data_clinical, y,
                         seed = 123, p_train = 0.7,
                         alpha_diff = 0.05) {

  X <- build_X_omics(merged_data)
  dc2 <- align_clinical(X, data_clinical)
  y2  <- align_y(X, data_clinical, y)

  y_bin <- ifelse(y2 == "no", 1, 0)

  set.seed(seed)
  train_idx <- createDataPartition(y2, p = p_train, list = FALSE)

  X_train_full <- X[train_idx, , drop = FALSE]
  X_test_full  <- X[-train_idx, , drop = FALSE]
  y_train <- y2[train_idx]
  y_test  <- y2[-train_idx]
  y_train_bin <- y_bin[train_idx]

  # remove zero-variance in TRAIN
  vars <- apply(X_train_full, 2, var)
  keep0 <- which(vars > 0 & !is.na(vars))
  X_train_full <- X_train_full[, keep0, drop = FALSE]
  X_test_full  <- X_test_full[,  keep0, drop = FALSE]

  # clinical (already dummies): all except ID
  clin_vars <- setdiff(colnames(dc2), "submitted_subject_id")
  clin_train <- dc2[train_idx, clin_vars, drop = FALSE]
  clin_test  <- dc2[-train_idx, clin_vars, drop = FALSE]

  # differential (TRAIN only) for ED__ and EXP__
  ed_cols  <- grep("^ED__",  colnames(X_train_full), value = TRUE)
  exp_cols <- grep("^EXP__", colnames(X_train_full), value = TRUE)

  p_ed  <- if (length(ed_cols)  > 0) sapply(ed_cols,  function(f) safe_p(X_train_full[[f]], y_train_bin)) else numeric(0)
  p_exp <- if (length(exp_cols) > 0) sapply(exp_cols, function(f) safe_p(X_train_full[[f]], y_train_bin)) else numeric(0)

  ed_keep  <- names(p_ed )[p_ed  < alpha_diff]
  exp_keep <- names(p_exp)[p_exp < alpha_diff]
  omic_keep <- c(ed_keep, exp_keep)

  # fallback to clinical-only if too few omics
  use_only_clinical <- FALSE
  if (length(omic_keep) < 2) {
    use_only_clinical <- TRUE
    X_train <- X_train_full[, 0, drop = FALSE]
    X_test  <- X_test_full[,  0, drop = FALSE]
  } else {
    X_train <- X_train_full[, omic_keep, drop = FALSE]
    X_test  <- X_test_full[,  omic_keep, drop = FALSE]
  }

  # LASSO (TRAIN only) if omics exist
  selected <- character(0)
  if (!use_only_clinical) {
    lasso_cv <- cv.glmnet(as.matrix(X_train), y_train_bin, alpha = 1, family = "binomial")
    coef_min <- coef(lasso_cv, s = "lambda.min")
    selected <- rownames(coef_min)[as.numeric(coef_min) != 0]
    selected <- setdiff(selected, "(Intercept)")

    if (length(selected) < 2) {
      use_only_clinical <- TRUE
      selected <- character(0)
      X_train_final <- clin_train
      X_test_final  <- clin_test
    } else {
      X_train_final <- cbind(X_train[, selected, drop = FALSE], clin_train)
      X_test_final  <- cbind(X_test[,  selected, drop = FALSE], clin_test)
    }
  } else {
    X_train_final <- clin_train
    X_test_final  <- clin_test
  }

  # models with CV in TRAIN
  ctrl <- trainControl(method = "cv", number = 10)

  train_df <- X_train_final
  train_df$y <- y_train

  set.seed(seed)
  model_rl <- train(y ~ ., data = train_df, method = "glmnet", family = "binomial", trControl = ctrl)

  set.seed(seed)
  model_rf <- train(y ~ ., data = train_df, method = "rf", trControl = ctrl)

  set.seed(seed)
  model_svm <- train(y ~ ., data = train_df, method = "svmLinear", trControl = ctrl)

  # evaluate on TEST
  pred_rl  <- predict(model_rl,  X_test_final)
  pred_rf  <- predict(model_rf,  X_test_final)
  pred_svm <- predict(model_svm, X_test_final)

  f1 <- c(
    RL  = unname(caret::confusionMatrix(pred_rl,  y_test, positive = "no")$byClass["F1"]),
    RF  = unname(caret::confusionMatrix(pred_rf,  y_test, positive = "no")$byClass["F1"]),
    SVM = unname(caret::confusionMatrix(pred_svm, y_test, positive = "no")$byClass["F1"])
  )

  list(
    seed = seed,
    F1 = f1,
    used_only_clinical = use_only_clinical,
    n_features_after_var0 = ncol(X_train_full),
    n_features_after_diff = ncol(X_train),
    n_features_lasso = ifelse(length(selected) == 0, 0, length(selected)),
    selected_omics = selected,
    y_test = y_test,
    pred_test = list(GLM = pred_rl, RF = pred_rf, SVM = pred_svm)
  )
}

##############################################################################
# 7) RUN MANY SEEDS -> all_runs
##############################################################################
stopifnot(exists("merged_list"), exists("data_clinical"), exists("selection"))
stopifnot(!is.null(selection$Breast.Nodal.CR.all))

seeds <- 1:50
alpha_diff <- 0.05

pb_seed <- progress_bar$new(
  format = "Seed [:bar] :current/:total (:percent) | elapsed :elapsed | eta :eta",
  total = length(seeds),
  clear = FALSE, width = 80
)
pb_seed$tick(0)

all_runs <- do.call(rbind, lapply(seeds, function(sd) {

  pb_seed$tick()

  pb_ds <- progress_bar$new(
    format = paste0("  Dataset [:bar] :current/:total (:percent) | elapsed :elapsed | eta :eta | seed ", sd),
    total = length(merged_list),
    clear = FALSE, width = 80
  )

  res_sd <- lapply(names(merged_list), function(ds) {
    pb_ds$tick()

    out <- run_pipeline(
      merged_data   = merged_list[[ds]],
      data_clinical = data_clinical,
      y             = selection$Breast.Nodal.CR.all,
      seed          = sd,
      p_train       = 0.7,
      alpha_diff    = alpha_diff
    )

    data.frame(
      seed = sd,
      Dataset = ds,
      Model = names(out$F1),
      F1 = as.numeric(out$F1),
      used_only_clinical = out$used_only_clinical,
      nfeat_var0 = out$n_features_after_var0,
      nfeat_diff = out$n_features_after_diff,
      nfeat_lasso = out$n_features_lasso,
      row.names = NULL
    )
  })

  do.call(rbind, res_sd)
}))

write.csv(all_runs, file.path("results", "metrics_all_runs_raw.csv"), row.names = FALSE)

##############################################################################
# 8) SUMMARY TABLE (mean + 95% CI)
##############################################################################
f1_summary <- all_runs %>%
  group_by(Dataset, Model) %>%
  summarise(
    n = n(),
    mean_F1 = mean(F1, na.rm = TRUE),
    sd_F1 = sd(F1, na.rm = TRUE),
    CI_low = quantile(F1, 0.025, na.rm = TRUE),
    CI_high = quantile(F1, 0.975, na.rm = TRUE),
    frac_only_clinical = mean(used_only_clinical, na.rm = TRUE),
    mean_nfeat_diff = mean(nfeat_diff, na.rm = TRUE),
    mean_nfeat_lasso = mean(nfeat_lasso, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(f1_summary, file.path("results", "metrics_summary_pretty.csv"), row.names = FALSE)

# Normalize model naming (RL -> GLM) for plots/tables
all_runs2 <- all_runs %>% mutate(Model = recode(Model, RL = "GLM"))

f1_summary2 <- all_runs2 %>%
  group_by(Dataset, Model) %>%
  summarise(
    n = sum(!is.na(F1)),
    mean_F1 = mean(F1, na.rm = TRUE),
    CI_low  = quantile(F1, 0.025, na.rm = TRUE),
    CI_high = quantile(F1, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(f1_summary2, file.path("results", "F1_summary_mean_CI95.csv"), row.names = FALSE)

##############################################################################
# 9) BARPLOT (mean ± 95% CI)  [base R, same as your code]
##############################################################################
f1_summary2$Model <- factor(f1_summary2$Model, levels = c("GLM", "RF", "SVM"))
f1_summary2$Dataset <- factor(f1_summary2$Dataset, levels = unique(all_runs2$Dataset))

mean_mat <- tapply(f1_summary2$mean_F1, list(f1_summary2$Model, f1_summary2$Dataset), c)
low_mat  <- tapply(f1_summary2$CI_low,  list(f1_summary2$Model, f1_summary2$Dataset), c)
high_mat <- tapply(f1_summary2$CI_high, list(f1_summary2$Model, f1_summary2$Dataset), c)

png(file.path("results", "F1_models_with_CI.png"), width = 1600, height = 900, res = 200)
op <- par(no.readonly = TRUE)
par(mar = c(10, 5, 4, 2) + 0.1, xpd = NA)

bp <- barplot(
  mean_mat,
  beside = TRUE,
  ylim = c(0, 1),
  ylab = "F1-score (mean ± 95% CI)",
  main = "Model performance across datasets",
  las = 2,
  legend.text = rownames(mean_mat),
  args.legend = list(x = "topright", bty = "n")
)

arrows(
  x0 = bp, y0 = low_mat,
  x1 = bp, y1 = high_mat,
  angle = 90, code = 3, length = 0.04
)

par(op)
dev.off()

##############################################################################
# 10) PICK REPRESENTATIVE SEED (median F1) + CONFUSION MATRIX (ED_EXP + GLM)
##############################################################################
seed_rep <- all_runs2 %>%
  filter(Dataset == "ED_EXP", Model == "GLM", !is.na(F1)) %>%
  mutate(dist = abs(F1 - median(F1, na.rm = TRUE))) %>%
  arrange(dist) %>%
  slice(1) %>%
  pull(seed)

out_one <- run_pipeline(
  merged_data   = merged_list$ED_EXP,
  data_clinical = data_clinical,
  y             = selection$Breast.Nodal.CR.all,
  seed          = seed_rep,
  p_train       = 0.7,
  alpha_diff    = 0.05
)

cm <- confusionMatrix(out_one$pred_test$GLM, out_one$y_test, positive = "no")
print(cm)
print(cm$byClass["F1"])

##############################################################################
# 11) FIGURE FINAL (A+B+C) (kept as your full block; only paths -> results/)
##############################################################################
# NOTE: This section preserves your logic; only file outputs go to ./results/

# ---- A) Paired ΔF1 with ED ----
comparisons <- tibble::tribble(
  ~Comparison,         ~base,      ~with_ed,
  "DNA \u2192 +ED",     "DNA",      "ED_DNA",
  "EXP \u2192 +ED",     "EXP",      "ED_EXP",
  "EXP+DNA \u2192 +ED", "EXP_DNA",  "ED_EXP_DNA"
)

pair_delta <- purrr::pmap_dfr(comparisons, function(Comparison, base, with_ed) {

  df_base <- all_runs2 %>% filter(Dataset == base) %>% select(seed, Model, F1_base = F1)
  df_ed   <- all_runs2 %>% filter(Dataset == with_ed) %>% select(seed, Model, F1_ed   = F1)

  full_join(df_base, df_ed, by = c("seed", "Model")) %>%
    mutate(Comparison = Comparison, deltaF1 = F1_ed - F1_base)
}) %>%
  filter(!is.na(deltaF1))

stats_A <- pair_delta %>%
  group_by(Model, Comparison) %>%
  summarise(
    n = n(),
    median_delta = median(deltaF1, na.rm = TRUE),
    p = suppressWarnings(wilcox.test(deltaF1, mu = 0, exact = FALSE)$p.value),
    pct_pos  = mean(deltaF1 > 0) * 100,
    pct_zero = mean(deltaF1 == 0) * 100,
    pct_neg  = mean(deltaF1 < 0) * 100,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf(
    "Wilcoxon p=%s\nMedian \u0394F1=%+.3f\n(+)%0.0f%% (=)%0.0f%% (-)%0.0f%% | n=%d",
    signif(p, 3), median_delta, pct_pos, pct_zero, pct_neg, n
  ))

ypos_A <- pair_delta %>%
  group_by(Model, Comparison) %>%
  summarise(y = max(deltaF1, na.rm = TRUE) + 0.03, .groups = "drop")

labels_A <- left_join(stats_A, ypos_A, by = c("Model", "Comparison")) %>%
  mutate(x = 1)

p_A <- ggplot(pair_delta, aes(x = 1, y = deltaF1)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.08, height = 0.01),
             alpha = 0.85, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.8) +
  facet_grid(Model ~ Comparison) +
  geom_text(data = labels_A, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.2) +
  scale_x_continuous(breaks = NULL) +
  labs(
    title = "A) Effect of adding RNA editing features (ED) on performance",
    subtitle = "Paired by seed. \u0394F1 = F1(with ED) \u2212 F1(without ED).",
    x = NULL,
    y = expression(Delta * "F1")
  ) +
  theme_classic(base_size = 12)

# ---- B) Feature importance (GLM) for representative seed (ED_EXP) ----
TOP_N <- 30
alpha_diff <- 0.05
cap_big_bar <- TRUE
cap_quantile <- 0.95

X  <- build_X_omics(merged_list$ED_EXP)
dc2 <- align_clinical(X, data_clinical)
y2  <- align_y(X, data_clinical, selection$Breast.Nodal.CR.all)

set.seed(seed_rep)
train_idx <- createDataPartition(y2, p = 0.7, list = FALSE)

X_train_full <- X[train_idx, , drop = FALSE]
y_train <- y2[train_idx]
y_bin   <- ifelse(y_train == "no", 1, 0)

vars <- apply(X_train_full, 2, var)
keep0 <- which(vars > 0 & !is.na(vars))
X_train_full <- X_train_full[, keep0, drop = FALSE]

clin_vars  <- setdiff(colnames(dc2), "submitted_subject_id")
clin_train <- dc2[train_idx, clin_vars, drop = FALSE]

ed_cols  <- grep("^ED__",  colnames(X_train_full), value = TRUE)
exp_cols <- grep("^EXP__", colnames(X_train_full), value = TRUE)

p_ed  <- if (length(ed_cols)  > 0) sapply(ed_cols,  function(f) safe_p(X_train_full[[f]], y_bin)) else numeric(0)
p_exp <- if (length(exp_cols) > 0) sapply(exp_cols, function(f) safe_p(X_train_full[[f]], y_bin)) else numeric(0)

omic_keep <- c(names(p_ed)[p_ed < alpha_diff], names(p_exp)[p_exp < alpha_diff])
if (length(omic_keep) < 2) stop("Representative seed left <2 omics after differential filtering.")

Xo <- X_train_full[, omic_keep, drop = FALSE]

lasso_cv <- cv.glmnet(as.matrix(Xo), y_bin, alpha = 1, family = "binomial")
coef_min <- coef(lasso_cv, s = "lambda.min")
selected <- setdiff(rownames(coef_min)[as.numeric(coef_min) != 0], "(Intercept)")
if (length(selected) < 1) stop("LASSO left 0 non-intercept features.")

X_final <- cbind(Xo[, selected, drop = FALSE], clin_train)

ctrl <- trainControl(method = "cv", number = 10)
df_train <- X_final
df_train$y <- y_train

set.seed(seed_rep)
model_glm <- train(
  y ~ ., data = df_train,
  method = "glmnet", family = "binomial",
  trControl = ctrl
)

coef_vec <- as.matrix(coef(model_glm$finalModel, s = model_glm$bestTune$lambda))

coef_df <- data.frame(
  Feature_raw = rownames(coef_vec),
  Beta = as.numeric(coef_vec[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(Feature_raw != "(Intercept)", Beta != 0) %>%
  mutate(
    Feature = trimws(gsub("`", "", Feature_raw)),
    Omics = case_when(
      grepl("^ED__",  Feature) ~ "RNA editing",
      grepl("^EXP__", Feature) ~ "Gene expression",
      TRUE ~ "Clinical"
    ),
    Feature_clean = gsub("^(ED__|EXP__)", "", Feature),
    AbsBeta = abs(Beta)
  ) %>%
  arrange(desc(AbsBeta))

if (cap_big_bar && nrow(coef_df) > 5) {
  cap <- as.numeric(quantile(coef_df$AbsBeta, cap_quantile, na.rm = TRUE))
  coef_df <- coef_df %>% mutate(AbsBeta_plot = pmin(AbsBeta, cap))
  y_lab <- paste0("|Coefficient| (capped at P", cap_quantile * 100, ")")
} else {
  coef_df <- coef_df %>% mutate(AbsBeta_plot = AbsBeta)
  y_lab <- "|Coefficient|"
}

coef_top <- coef_df %>% slice_head(n = TOP_N)
coef_top$Omics <- factor(coef_top$Omics, levels = c("Clinical", "Gene expression", "RNA editing"))

p_B <- ggplot(coef_top, aes(x = reorder(Feature_clean, AbsBeta_plot),
                            y = AbsBeta_plot, fill = Omics)) +
  geom_col(width = 0.85) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Clinical" = "#4D4D4D",
    "Gene expression" = "#1F77B4",
    "RNA editing" = "#D62728"
  )) +
  labs(
    title = paste0("B) Top GLM features (seed ", seed_rep, ")"),
    subtitle = "Importance = |glmnet coefficient| (non-zero at best lambda)",
    x = NULL, y = y_lab, fill = "Feature type"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"))

# ---- C) Confusion matrix plot (TEST) ----
label_map <- c("yes" = "Responder", "no" = "Non-responder")
level_order <- c("Non-responder", "Responder")

cm_df <- as.data.frame(cm$table) %>%
  rename(Predicted_raw = Prediction, Observed_raw = Reference, n = Freq) %>%
  mutate(
    Predicted = recode(as.character(Predicted_raw), !!!label_map),
    Observed  = recode(as.character(Observed_raw),  !!!label_map),
    Predicted = factor(Predicted, levels = level_order),
    Observed  = factor(Observed,  levels = level_order)
  )

acc  <- unname(cm$overall["Accuracy"])
sens <- unname(cm$byClass["Sensitivity"])
spec <- unname(cm$byClass["Specificity"])
f1no <- unname(cm$byClass["F1"])

p_C <- ggplot(cm_df, aes(x = Observed, y = Predicted)) +
  geom_tile(fill = "white", color = "black", linewidth = 0.8) +
  geom_text(aes(label = n), size = 10, fontface = "bold") +
  coord_equal() +
  labs(
    title = "C) Confusion matrix (test set)",
    subtitle = sprintf("Positive class = Non-responder ('no') | Acc=%.2f | Sens=%.2f | Spec=%.2f | F1(no)=%.2f",
                       acc, sens, spec, f1no),
    x = "Observed",
    y = "Predicted"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

theme_global <- theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text    = element_text(size = 12, face = "bold"),
    axis.title    = element_text(size = 13),
    axis.text     = element_text(size = 11),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 11),
    plot.margin   = ggplot2::margin(t = 8, r = 10, b = 8, l = 10, unit = "pt")
  )

final <- (p_A | p_B) +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 18, face = "bold"))
  ) &
  theme_global

print(final)

ggsave(file.path("results", "Figure_final_A_B_C_seedRep.png"),
       final, width = 14, height = 8, units = "in", dpi = 400, limitsize = FALSE)
ggsave(file.path("results", "Figure_final_A_B_C_seedRep.pdf"),
       final, width = 14, height = 8, units = "in", limitsize = FALSE)

cat("\nSaved:\n- results/Figure_final_A_B_C_seedRep.png\n- results/Figure_final_A_B_C_seedRep.pdf\n")
cat("Seed used =", seed_rep, "\n")
cat("F1(no) in TEST =", f1no, "\n")

##############################################################################
# 12) TABLE: variables in Figure B (Top coefficients)
##############################################################################
coef_df2 <- data.frame(
  Feature_raw = rownames(coef_vec),
  Beta = as.numeric(coef_vec[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(Feature_raw != "(Intercept)", Beta != 0) %>%
  mutate(
    Feature = gsub("`", "", Feature_raw),
    Omics = case_when(
      grepl("^ED__",  Feature) ~ "RNA editing",
      grepl("^EXP__", Feature) ~ "Gene expression",
      TRUE ~ "Clinical"
    ),
    Feature = gsub("^(ED__|EXP__)", "", Feature),
    OR = exp(Beta),
    Direction = ifelse(Beta > 0, "Higher non-responder risk", "Lower non-responder risk"),
    AbsBeta = abs(Beta)
  ) %>%
  arrange(desc(AbsBeta)) %>%
  mutate(
    Rank = row_number(),
    Beta = formatC(Beta, format = "f", digits = 4),
    OR   = formatC(OR,   format = "f", digits = 3)
  ) %>%
  select(Rank, Feature, Omics, Beta, OR, Direction)

write.csv(coef_df2, file.path("results", "Fig3B_top_features.csv"), row.names = FALSE)
write.table(coef_df2, file.path("results", "Fig3B_top_features.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

##############################################################################
# 13) EXTRA METRICS (MCC, PR-AUC) — ED_EXP + GLM — positive="no"
##############################################################################
mcc_from_cm <- function(cm_table, positive = "no") {
  classes <- colnames(cm_table)
  stopifnot(all(classes %in% rownames(cm_table)))
  stopifnot(positive %in% classes)
  neg <- setdiff(classes, positive)
  stopifnot(length(neg) == 1)
  neg <- neg[[1]]

  TP <- cm_table[positive, positive]
  TN <- cm_table[neg, neg]
  FP <- cm_table[positive, neg]
  FN <- cm_table[neg, positive]

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (denom == 0) return(NA_real_)
  as.numeric((TP * TN - FP * FN) / denom)
}

set.seed(seed_rep)
ctrl_prob <- trainControl(method = "cv", number = 10, classProbs = TRUE, savePredictions = "final")

alpha_diff <- 0.05
p_train <- 0.7

X  <- build_X_omics(merged_list$ED_EXP)
dc2 <- align_clinical(X, data_clinical)
y2  <- align_y(X, data_clinical, selection$Breast.Nodal.CR.all)

set.seed(seed_rep)
train_idx <- createDataPartition(y2, p = p_train, list = FALSE)

X_train_full <- X[train_idx, , drop = FALSE]
X_test_full  <- X[-train_idx, , drop = FALSE]
y_train <- y2[train_idx]
y_test  <- y2[-train_idx]
y_train_bin <- ifelse(y_train == "no", 1, 0)

vars <- apply(X_train_full, 2, var)
keep0 <- which(vars > 0 & !is.na(vars))
X_train_full <- X_train_full[, keep0, drop = FALSE]
X_test_full  <- X_test_full[,  keep0, drop = FALSE]

clin_vars  <- setdiff(colnames(dc2), "submitted_subject_id")
clin_train <- dc2[train_idx, clin_vars, drop = FALSE]
clin_test  <- dc2[-train_idx, clin_vars, drop = FALSE]

ed_cols  <- grep("^ED__",  colnames(X_train_full), value = TRUE)
exp_cols <- grep("^EXP__", colnames(X_train_full), value = TRUE)

p_ed  <- if (length(ed_cols)  > 0) sapply(ed_cols,  function(f) safe_p(X_train_full[[f]], y_train_bin)) else numeric(0)
p_exp <- if (length(exp_cols) > 0) sapply(exp_cols, function(f) safe_p(X_train_full[[f]], y_train_bin)) else numeric(0)

omic_keep <- c(names(p_ed)[p_ed < alpha_diff], names(p_exp)[p_exp < alpha_diff])

use_only_clinical <- FALSE
if (length(omic_keep) < 2) {
  use_only_clinical <- TRUE
  X_train_final <- clin_train
  X_test_final  <- clin_test
} else {
  Xo_train <- X_train_full[, omic_keep, drop = FALSE]
  Xo_test  <- X_test_full[,  omic_keep, drop = FALSE]

  lasso_cv <- cv.glmnet(as.matrix(Xo_train), y_train_bin, alpha = 1, family = "binomial")
  coef_min <- coef(lasso_cv, s = "lambda.min")
  selected <- rownames(coef_min)[as.numeric(coef_min) != 0]
  selected <- setdiff(selected, "(Intercept)")

  if (length(selected) < 2) {
    use_only_clinical <- TRUE
    X_train_final <- clin_train
    X_test_final  <- clin_test
  } else {
    X_train_final <- cbind(Xo_train[, selected, drop = FALSE], clin_train)
    X_test_final  <- cbind(Xo_test[,  selected, drop = FALSE], clin_test)
  }
}

train_df <- X_train_final
train_df$y <- y_train

set.seed(seed_rep)
model_glm_prob <- train(
  y ~ ., data = train_df,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl_prob
)

pred_class <- predict(model_glm_prob, X_test_final)
pred_prob  <- predict(model_glm_prob, X_test_final, type = "prob")[, "no"]

cm2 <- confusionMatrix(pred_class, y_test, positive = "no")
mcc <- mcc_from_cm(cm2$table, positive = "no")

pr <- PRROC::pr.curve(
  scores.class0 = pred_prob[y_test == "no"],
  scores.class1 = pred_prob[y_test == "yes"],
  curve = FALSE
)

cat("\n================ EXTRA METRICS (ED_EXP + GLM) ================\n")
cat("Seed:", seed_rep, "\n")
cat("Clinical-only fallback?:", use_only_clinical, "\n\n")
print(cm2)
cat("\nMCC =", round(mcc, 4), "\n")
cat("PR-AUC (positive='no') =", round(pr$auc.integral, 4), "\n")
cat("F1 (positive='no') =", round(unname(cm2$byClass["F1"]), 4), "\n")
cat("=============================================================\n")

##############################################################################
# 14) (Optional) Metrics per seed (ED_EXP + GLM): F1, MCC, PR-AUC
##############################################################################
run_one_glm_with_prob <- function(merged_data, data_clinical, y,
                                  seed = 123, p_train = 0.7, alpha_diff = 0.05) {

  X  <- build_X_omics(merged_data)
  dc2 <- align_clinical(X, data_clinical)
  y2  <- align_y(X, data_clinical, y)
  y_bin <- ifelse(y2 == "no", 1, 0)

  set.seed(seed)
  train_idx <- createDataPartition(y2, p = p_train, list = FALSE)

  X_train_full <- X[train_idx, , drop = FALSE]
  X_test_full  <- X[-train_idx, , drop = FALSE]
  y_train <- y2[train_idx]
  y_test  <- y2[-train_idx]
  y_train_bin <- y_bin[train_idx]

  vars <- apply(X_train_full, 2, var)
  keep0 <- which(vars > 0 & !is.na(vars))
  X_train_full <- X_train_full[, keep0, drop = FALSE]
  X_test_full  <- X_test_full[,  keep0, drop = FALSE]

  clin_vars  <- setdiff(colnames(dc2), "submitted_subject_id")
  clin_train <- dc2[train_idx, clin_vars, drop = FALSE]
  clin_test  <- dc2[-train_idx, clin_vars, drop = FALSE]

  ed_cols  <- grep("^ED__",  colnames(X_train_full), value = TRUE)
  exp_cols <- grep("^EXP__", colnames(X_train_full), value = TRUE)

  p_ed  <- if (length(ed_cols)  > 0) sapply(ed_cols,  function(f) safe_p(X_train_full[[f]], y_train_bin)) else numeric(0)
  p_exp <- if (length(exp_cols) > 0) sapply(exp_cols, function(f) safe_p(X_train_full[[f]], y_train_bin)) else numeric(0)

  omic_keep <- c(names(p_ed)[p_ed < alpha_diff], names(p_exp)[p_exp < alpha_diff])

  if (length(omic_keep) < 2) {
    X_train_final <- clin_train
    X_test_final  <- clin_test
  } else {
    Xo_train <- X_train_full[, omic_keep, drop = FALSE]
    Xo_test  <- X_test_full[,  omic_keep, drop = FALSE]

    lasso_cv <- cv.glmnet(as.matrix(Xo_train), y_train_bin, alpha = 1, family = "binomial")
    coef_min <- coef(lasso_cv, s = "lambda.min")
    selected <- rownames(coef_min)[as.numeric(coef_min) != 0]
    selected <- setdiff(selected, "(Intercept)")

    if (length(selected) < 2) {
      X_train_final <- clin_train
      X_test_final  <- clin_test
    } else {
      X_train_final <- cbind(Xo_train[, selected, drop = FALSE], clin_train)
      X_test_final  <- cbind(Xo_test[,  selected, drop = FALSE], clin_test)
    }
  }

  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, savePredictions = "final")

  train_df <- X_train_final
  train_df$y <- y_train

  set.seed(seed)
  model_glm <- train(
    y ~ ., data = train_df,
    method = "glmnet",
    family = "binomial",
    trControl = ctrl
  )

  pred_class <- predict(model_glm, X_test_final)
  pred_prob  <- predict(model_glm, X_test_final, type = "prob")[, "no"]

  list(
    seed = seed,
    y_test = y_test,
    pred_test = list(GLM = pred_class),
    prob_test = list(GLM = pred_prob)
  )
}

metrics_list <- lapply(seeds, function(sd) {

  out <- run_one_glm_with_prob(
    merged_data   = merged_list$ED_EXP,
    data_clinical = data_clinical,
    y             = selection$Breast.Nodal.CR.all,
    seed          = sd,
    p_train       = 0.7,
    alpha_diff    = 0.05
  )

  cm <- confusionMatrix(out$pred_test$GLM, out$y_test, positive = "no")
  tab <- cm$table

  TP <- tab["no",  "no"]
  TN <- tab["yes", "yes"]
  FP <- tab["no",  "yes"]
  FN <- tab["yes", "no"]

  mcc <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  prob <- out$prob_test$GLM
  y_test <- out$y_test
  pr <- PRROC::pr.curve(
    scores.class0 = prob[y_test == "no"],
    scores.class1 = prob[y_test == "yes"],
    curve = FALSE
  )

  data.frame(
    seed = sd,
    F1 = unname(cm$byClass["F1"]),
    MCC = as.numeric(mcc),
    PR_AUC = as.numeric(pr$auc.integral),
    stringsAsFactors = FALSE
  )
})

metrics_all <- dplyr::bind_rows(metrics_list)
write.csv(metrics_all, file.path("results", "ED_EXP_GLM_metrics_per_seed.csv"), row.names = FALSE)

metrics_summary <- metrics_all %>%
  summarise(
    n = n(),
    mean_F1 = mean(F1, na.rm = TRUE),
    sd_F1 = sd(F1, na.rm = TRUE),
    mean_MCC = mean(MCC, na.rm = TRUE),
    sd_MCC = sd(MCC, na.rm = TRUE),
    mean_PR_AUC = mean(PR_AUC, na.rm = TRUE),
    sd_PR_AUC = sd(PR_AUC, na.rm = TRUE)
  )

print(metrics_summary)

# Optional plots (same as your idea)
metrics_long <- metrics_all %>%
  select(seed, MCC, PR_AUC) %>%
  pivot_longer(cols = c(MCC, PR_AUC), names_to = "metric", values_to = "value")

sum_long <- metrics_long %>%
  group_by(metric) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), .groups = "drop")

p_metrics <- ggplot(metrics_long, aes(x = seed, y = value)) +
  geom_line() +
  geom_point(size = 1.6) +
  geom_hline(data = sum_long, aes(yintercept = mean), inherit.aes = FALSE, linetype = 2) +
  geom_hline(data = sum_long, aes(yintercept = mean + sd), inherit.aes = FALSE, linetype = 3) +
  geom_hline(data = sum_long, aes(yintercept = mean - sd), inherit.aes = FALSE, linetype = 3) +
  facet_grid(metric ~ ., scales = "free_y") +
  labs(x = "Seed", y = "Value", title = "MCC and PR-AUC per seed in ED_EXP/GLM") +
  theme_bw()

ggsave(file.path("results", "ED_EXP_GLM_MCC_PR_AUC_by_seed.png"),
       p_metrics, width = 10, height = 6, dpi = 250)
