# Required Libraries
library(dplyr)
library(tidyr)
library(biomaRt)

# Directories
base_path <- "data/input/20240605_beauty_salmon_results/"
clinical_path <- "data/clinical/merge_final_104.tsv"
phs_path <- "data/input/phs001050.v1.p1.txt"
beauty_path <- "data/input/BEAUTY_SraRunTable.txt"
pivoted_output <- "output/all_data_pivoted_104_DEG.tsv"
genes_output <- "output/all_data_pivoted_with_genes.tsv"

# IDs for analysis
ids <- 3069822:3070020
data_list <- list()

# Read and process files
for (id in ids) {
  file <- file.path(base_path, paste0("SRR", id, "/quant.sf"))
  if (file.exists(file)) {
    df <- read.delim(file, select = c("Name", "NumReads"))
    df$Sample <- paste0("SRR", id)
    data_list[[paste0("SRR", id)]] <- df
  } else {
    warning(paste("File not found:", file))
  }
}

# Combine and pivot data
all_data <- bind_rows(data_list)
all_data_pivoted <- all_data %>%
  pivot_wider(names_from = Sample, values_from = NumReads)

# Clean up column names
colnames(all_data_pivoted) <- gsub("^.*_", "", colnames(all_data_pivoted))

# Load and process additional data
merge_final_104 <- read.csv(clinical_path, sep = "\t")
phs001050 <- read.delim(phs_path) %>%
  filter(SUBJECT_ID %in% merge_final_104$submitted_subject_id, sequence_sample_use == "RNA-Seq")
BEAUTY_SraRunTable <- read.csv(beauty_path) %>%
  filter(Assay.Type == "RNA-Seq")

# Filter and merge relevant data
merged_data <- phs001050 %>%
  inner_join(BEAUTY_SraRunTable, by = c("SUBJECT_ID" = "submitted_subject_id"))
all_data_pivoted <- all_data_pivoted %>%
  select(Name, all_of(merged_data$Run)) %>%
  as.data.frame()

# Save intermediate result
write.table(all_data_pivoted, pivoted_output, row.names = FALSE, sep = "\t")

# Replace Ensembl IDs with HUGO symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
hgnc_symbols <- getBM(
  attributes = c("ensembl_transcript_id_version", "external_gene_name"),
  filters = "ensembl_transcript_id_version",
  values = rownames(all_data_pivoted),
  mart = ensembl
)

# Map gene names
all_data_pivoted <- all_data_pivoted %>%
  rownames_to_column("Name") %>%
  left_join(hgnc_symbols, by = c("Name" = "ensembl_transcript_id_version")) %>%
  relocate(external_gene_name, .before = everything()) %>%
  select(-Name)

# Rename primary column
colnames(all_data_pivoted)[1] <- "Gene"

# Save final result
write.table(all_data_pivoted, genes_output, row.names = FALSE, sep = "\t")

# Finalization message
cat("Data processing is complete. Files have been saved in 'output/'.\n")



################## VOLCANO PLOT #################################################



# Set the working directory
setwd("path_to_working_directory/")

# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("EnhancedVolcano", force = TRUE)
library(DESeq2)
library(dplyr)
library(biomaRt)
library(EnhancedVolcano)

# Read and merge data
merge <- read.csv("path_to_merge_final_104.tsv", sep = "\t")

merged_data_1 <- merged_data %>%
  inner_join(merge, by = c("SUBJECT_ID" = "submitted_subject_id"))

# Select relevant columns and rename
merge <- dplyr::select(merged_data_1, "Run.x", "Breast.Nodal.CR.all")
merge$condition <- merge$Breast.Nodal.CR.all
merge$Breast.Nodal.CR.all <- NULL

# Prepare coldata for DESeq2
coldata <- as.data.frame(merge)
rownames(coldata) <- coldata$Run

# Prepare countdata for DESeq2
countdata <- as.data.frame(all_data_pivoted)
countdata <- countdata[, coldata$Run.x]
countdata_rounded <- round(countdata, digits = 0)
countdata_integer <- as.matrix(countdata_rounded)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata_integer,
                              colData = coldata,
                              design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)

# Retrieve results
res <- results(dds)

# Annotate results with gene symbols using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
hgnc_symbols <- getBM(attributes = c("ensembl_transcript_id_version", "hgnc_symbol"),
                      filters = "ensembl_transcript_id_version",
                      values = rownames(res),
                      mart = ensembl)
res$hugosymbol <- hgnc_symbols$hgnc_symbol[match(rownames(res), hgnc_symbols$ensembl_transcript_id_version)]


# MA-plot
plotMA(res, main = "DESeq2", ylim = c(-30, 30))

# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

# Filter significant genes
significant_genes <- rownames(res)[res$padj < 0.05 & abs(res$log2FoldChange) > log2(2.5)]
filtered_data <- all_data_pivoted[rownames(all_data_pivoted) %in% significant_genes, ]


# Volcano plot with enhanced labeling
nombres_genes_modificados <- ifelse(rownames(res) %in% significant_genes, res$hugosymbol, rownames(res))
EnhancedVolcano(res,
                lab = nombres_genes_modificados,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 2.5,
                pointSize = 0.5,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'top',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                max.overlaps = 10)

# Print confirmation
print("Analysis complete!")

