############################## tumor DNA analysis ################################################

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)

# Define working directory and path to MAF files
setwd("path_to_working_directory/")
folder_path <- "path_to_maf_files/"

# List all MAF files in the folder
maf_files <- list.files(path = folder_path, pattern = "\\.maf$", full.names = TRUE)

# Function to read MAF files and add sample names
read_maf_file <- function(file) {
  maf_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  sample_name <- gsub("\\.filtered\\.vep\\.vcf2maf\\.maf$", "", basename(file))
  maf_data$Sample <- sample_name
  return(maf_data)
}

# Read and combine all MAF files
maf_list <- lapply(maf_files, read_maf_file)
combined_maf <- do.call(rbind, maf_list)

# Filter MAF data to include only samples in `merge_filtered$Run`
filtered_maf <- combined_maf %>%
  filter(Tumor_Sample_Barcode %in% merge_filtered$Run)

# Define DDR genes
genes_ddr <- c("APLF", "APTX", "ASCC3", "DNTT", "LIG1", "LIG3", "LIG4", "MRE11",
               "NBN", "NHEJ1", "PARG", "PARP1", "PARP3", "PARPBP", "PNKP", "POLB",
               # (Include all other DDR genes here)
               "XPC", "ZSWIM7", "PTEN", "TDP2", "ENDOV", "SPRTN", "RNF4", 
               "SMARCA4", "IDH1", "SOX4", "WEE1", "RAD9B", "AEN", "PLK3")

# Filter MAF for DDR genes
filtered_genes <- filtered_maf %>%
  filter(Hugo_Symbol %in% genes_ddr, FILTER == "PASS")

# Add a unique identifier for each mutation
filtered_genes <- filtered_genes %>%
  mutate(tDNA_ID = paste("tDNA", Chromosome, Start_Position, End_Position, Strand, Hugo_Symbol, sep = "_"))

# Remove duplicate entries and prepare for pivoting
maf_pass_g_list <- filtered_genes %>%
  distinct(tDNA_ID, Tumor_Sample_Barcode, .keep_all = TRUE)

# Pivot data: Create a matrix with mutations
maf_pivot <- maf_pass_g_list %>%
  select(tDNA_ID, Tumor_Sample_Barcode) %>%
  mutate(has_mutation = 1) %>%
  pivot_wider(id_cols = tDNA_ID, names_from = Tumor_Sample_Barcode, values_from = has_mutation,
              values_fill = list(has_mutation = 0))

# Load gene consensus file
geneconsensus <- read.delim("path_to_geneconsensus_file.tsv")

# Filter maf_pivot for genes in the consensus list and DDR genes
maf_pivot <- maf_pivot %>%
  filter(str_extract(tDNA_ID, "[^_]+$") %in% geneconsensus$Gene.Symbol)

maf_pivot_ddr <- maf_pivot %>%
  filter(str_extract(tDNA_ID, "[^_]+$") %in% genes_ddr)

# Map column names to submitted subject IDs
mapping <- merge_filtered %>%
  select(Run, submitted_subject_id) %>%
  distinct()

current_names <- colnames(maf_pivot)
current_names_tumor <- current_names[-1]

new_names <- sapply(current_names_tumor, function(x) {
  if (x %in% mapping$Run) {
    mapping$submitted_subject_id[match(x, mapping$Run)]
  } else {
    x
  }
})

colnames(maf_pivot) <- c("tDNA_ID", new_names)
colnames(maf_pivot_ddr) <- c("tDNA_ID", new_names)

# Preview results
head(maf_pivot)
head(maf_pivot_ddr)


####################### germline DNA analysis #######################

# Load required libraries
library(dplyr)
library(tidyr)
library(maftools)
library(RColorBrewer)

# Define working directory and path to MAF files
folder_path <- "path_to_normal_sample_vep_mafs/"
maf_files <- list.files(path = folder_path, pattern = "\\.maf$", full.names = TRUE)

# Function to read MAF files and add sample names
read_maf_file <- function(file) {
  maf_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  sample_name <- gsub("\\.filtered\\.vep\\.vcf2maf\\.maf$", "", basename(file))
  maf_data$Sample <- sample_name
  return(maf_data)
}

# Read and combine all MAF files
maf_list <- lapply(maf_files, read_maf_file)
combined_maf_g <- do.call(rbind, maf_list)

# Clean the Sample column
combined_maf_g$Sample <- gsub("\\.maf$", "", combined_maf_g$Sample)

# Load DEG and clinical data
deg <- read.delim("path_to_deg_file/v1_all_data_pivoted_filtered_104_subject_id.txt")
rownames(deg) <- deg$Name
deg$Name <- NULL
Run1 <- colnames(deg)

clinical_data_forR <- read.delim("path_to_clinical_data/clinical_data_forR.tsv")
BEAUTY_SraRunTable <- read.csv("path_to_BEAUTY_SraRunTable.csv")
phs001050.v1.p1 <- read.delim("path_to_phs001050_file.txt")

# Merge clinical data
merge <- merge(BEAUTY_SraRunTable, clinical_data_forR, by.x = "submitted_subject_id", by.y = "Subject_ID", all = TRUE)
merge_filtered <- subset(merge, submitted_subject_id %in% Run1)
merge_filtered <- filter(merge_filtered, analyte_type == "DNA" & body_site == "Blood")

# Filter combined MAF for relevant samples
combined_maf_g <- subset(combined_maf_g, Sample %in% merge_filtered$Run)
maf_pass_g <- filter(combined_maf_g, FILTER == "PASS")

# Define the gene list
gene_list <- c("ATM", "BAPI", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "MSH2", "MSH6",
               "MUTYH", "DICER1", "PALB2", "RUNX1", "SDHAF2", "SDHB", "SDHC",
               # Add other genes as required
               "RECQL4", "TP53")

# Filter for genes in the list and add unique gDNA_ID
maf_pass_g_list <- maf_pass_g %>%
  filter(Hugo_Symbol %in% gene_list) %>%
  mutate(gDNA_ID = paste("gDNA", Chromosome, Start_Position, End_Position, Strand, Hugo_Symbol, sep = "_")) %>%
  distinct(gDNA_ID, Tumor_Sample_Barcode, .keep_all = TRUE)

# Pivot data
maf_pivot <- maf_pass_g_list %>%
  select(gDNA_ID, Tumor_Sample_Barcode) %>%
  mutate(has_mutation = 1) %>%
  pivot_wider(id_cols = gDNA_ID, names_from = Tumor_Sample_Barcode, values_from = has_mutation,
              values_fill = list(has_mutation = 0))

write.table(maf_pivot, "maf_pivoted_gDNA_allgenes104.tsv", row.names = FALSE)

# Filter for missense variants
maf_pass_g_missense <- filter(maf_pass_g_list, Consequence == "missense_variant" & CLIN_SIG != "benign") %>%
  distinct(gDNA_ID, Tumor_Sample_Barcode, .keep_all = TRUE)

maf_pivot_missense <- maf_pass_g_missense %>%
  select(gDNA_ID, Tumor_Sample_Barcode) %>%
  mutate(has_mutation = 1) %>%
  pivot_wider(id_cols = gDNA_ID, names_from = Tumor_Sample_Barcode, values_from = has_mutation,
              values_fill = list(has_mutation = 0))

write.table(maf_pivot_missense, "gDNA_missense_pivoted104.tsv", row.names = FALSE)

# Rename columns in maf_pivot
mapping <- merge_filtered %>%
  select(Run, submitted_subject_id) %>%
  distinct()

current_names <- colnames(maf_pivot)
current_names_tumor <- current_names[-1]

new_names <- sapply(current_names_tumor, function(x) {
  if (x %in% mapping$Run) {
    mapping$submitted_subject_id[match(x, mapping$Run)]
  } else {
    x
  }
})

colnames(maf_pivot) <- c("gDNA_ID", new_names)
colnames(maf_pivot_missense) <- c("gDNA_ID", new_names)

# Visualize results
head(maf_pivot)
head(maf_pivot_missense)

# Filter and create MAF object
required_columns <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                      "Variant_Classification", "Variant_Type", "Reference_Allele",
                      "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")

maf_data <- maf_pass_g_list[, required_columns] %>%
  na.omit()

maf_object <- read.maf(maf = maf_data)

# Add clinical data to the MAF object
merge_filtered$Run <- trimws(merge_filtered$Run)
maf_object_filtered <- read.maf(maf = maf_data, clinicalData = merge_filtered)

# Generate oncoplot
oncoplot(maf = maf_object_filtered,
         top = 10,
         showTumorSampleBarcodes = FALSE)

# Create coOncoplot for groups
maf_no <- subsetMaf(maf = maf_object_filtered, tsb = maf_object_filtered@clinical.data$Tumor_Sample_Barcode[maf_object_filtered@clinical.data[["Breast.Nodal.CR.all"]] == "no"])
maf_yes <- subsetMaf(maf = maf_object_filtered, tsb = maf_object_filtered@clinical.data$Tumor_Sample_Barcode[maf_object_filtered@clinical.data[["Breast.Nodal.CR.all"]] == "yes"])

coOncoplot(
  m1 = maf_no,
  m2 = maf_yes,
  m1Name = "No",
  m2Name = "Yes",
  clinicalFeatures1 = "Clinical.Molecular.Subtype",
  clinicalFeatures2 = "Clinical.Molecular.Subtype"
)

