# Load required libraries
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)

# Define working directory and input file
setwd("path_to_input_directory/")
reditools <- read.delim("path_to_input_file/BEAUTY107_master_table_biallelic_minimal.tsv")
reditools <- as.data.table(reditools)

# Select relevant columns
reditools <- reditools %>%
  select(c(1, 2, 3, 97:ncol(reditools))) %>%
  filter(!if_all(4:ncol(reditools), is.na)) %>%
  filter(Reference == "A" | Reference == "T")

# Identify base count columns
base_count_columns <- grep("BaseCount.A.C.G.T.", names(reditools), value = TRUE)

# Create a new dataset with selected columns
new_dataset <- reditools[, .(Region, Position, Reference)]

# Add columns for reference and alternate base counts
for (sample_col in base_count_columns) {
  sample_name <- sub("BaseCount.A.C.G.T._", "", sample_col)
  
  new_dataset[, paste0("BaseCount_ref_", sample_name) := sapply(seq_len(nrow(reditools)), function(i) {
    counts <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", reditools[i, get(sample_col)]), ",")))
    if (length(counts) == 4) {
      ifelse(reditools$Reference[i] == "A", counts[1], ifelse(reditools$Reference[i] == "T", counts[4], NA))
    } else {
      NA
    }
  })]
  
  new_dataset[, paste0("BaseCount_alt_", sample_name) := sapply(seq_len(nrow(reditools)), function(i) {
    counts <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", reditools[i, get(sample_col)]), ",")))
    if (length(counts) == 4) {
      ifelse(reditools$Reference[i] == "A", counts[3], ifelse(reditools$Reference[i] == "T", counts[2], NA))
    } else {
      NA
    }
  })]
}

# Add Uploaded_variation column
reditools[, Uploaded_variation := paste0(Region, "_", Position, ifelse(Reference == "A", "_A/G", ifelse(Reference == "T", "_T/C", "")))]

# Filter for relevant SRR columns
srr_cols <- grep("SRR", names(reditools), value = TRUE)
srr_cols <- c("Uploaded_variation", srr_cols)
srr_data <- reditools[, ..srr_cols]

# Rename SRR columns to keep only the SRR number
new_names <- gsub(".*_(SRR\\d+).*", "\\1", srr_cols)
setnames(srr_data, old = srr_cols[-1], new = new_names[-1])


# Load VEP (variant effect predictor) results and filter
vep <- read.delim("path_to_output_file/vep.txt")
vep_filter <- filter(vep, Existing_variation == "-")
ID <- vep_filter$Uploaded_variation
srr_filter$Uploaded_variation <- gsub("chr", "", srr_filter$Uploaded_variation)
srr_filter_final <- srr_filter[srr_filter$Uploaded_variation %in% ID, , drop = FALSE]
write.table(srr_filter_final, "reditool_final.csv", row.names = FALSE)

# Apply statistical model functions
source("REDIT_LLR.R")
out.df <- data.frame()
for (i in 1:nrow(srr_data)) {
  mat.ed <- create.matrix(srr_data[i, 2:209])
  mat.ed <- mat.ed[, colSums(is.na(mat.ed)) == 0]
  group <- ifelse(colnames(mat.ed) %in% non_responder, "NR", "R")
  group <- as.factor(group)
  mat.ed <- matrix(as.numeric(mat.ed), ncol = length(colnames(mat.ed)))
  test.tmp <- REDIT_LLR(data = mat.ed, groups = group)
  out.df <- rbind(out.df, data.frame(
    srr_data[i, ],
    pvalue = test.tmp$p.value,
    mle_group_resistance_alpha = test.tmp$mle.for.group.NR[[1]],
    mle_group_resistance_beta = test.tmp$mle.for.group.NR[[2]],
    mle_group_sensitive_alpha = test.tmp$mle.for.group.R[[1]],
    mle_group_sensitive_beta = test.tmp$mle.for.group.R[[2]],
    mle_null_model = test.tmp$mle.for.null.model,
    log_l_resistance = test.tmp$log.likelihood.for.group.NR,
    log_l_sensitive = test.tmp$log.likelihood.for.group.R,
    log_l_null = test.tmp$log.likelihood.for.null
  ))
}

# Calculate and annotate significant results
data <- out.df[!duplicated(out.df$Uploaded_variation), ]
data$EL_NR <- data$mle_group_resistance_alpha / (data$mle_group_resistance_alpha + data$mle_group_resistance_beta)
data$EL_R <- data$mle_group_sensitive_alpha / (data$mle_group_sensitive_alpha + data$mle_group_sensitive_beta)
data$delta <- data$EL_R - data$EL_NR
data$class <- ifelse(data$delta >= 0, "responder", "non_responder")
data$fdr <- p.adjust(data$pvalue, method = "BH")
data_significant <- filter(data, fdr < 0.01 & abs(delta) > 0.05)

# Generate volcano plot
ggplot(data, aes(x = delta, y = -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr < 0.01 & abs(delta) > 0.05, "significant", "not_significant")), alpha = 0.6, size = 2) +
  labs(title = "Volcano Plot", x = "Delta (EL_NR - EL_R)", y = "-log10(FDR)") +
  theme_minimal() +
  scale_color_manual(values = c("significant" = "#009E73", "not_significant" = "gray")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", size = 0.8) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", size = 0.8) +
  geom_vline(xintercept = -0.05, linetype = "dashed", color = "black", size = 0.8) +
  theme(legend.title = element_blank())
