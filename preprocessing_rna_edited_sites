# Load required libraries
library(tidyr)
library(dplyr)
library(data.table)

# Define working directory and input file
setwd("path_to_input_directory/")
reditools <- read.delim("path_to_input_file/BEAUTY107_master_table_biallelic_minimal.tsv")
reditools <- as.data.table(reditools)

# Identify columns containing base counts
base_count_columns <- grep("BaseCount.A.C.G.T", names(reditools), value = TRUE)

# Process base count columns
for (col in base_count_columns) {
  reditools[, paste0(col, "_Value") := sapply(get(col), function(counts) {
    # Convert BaseCount[A,C,G,T] values to numeric
    counts <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", counts), ",")))
    # Handle values based on the Reference column
    if (length(counts) == 4) {
      if (Reference[1] == "A") {
        return(counts[1] / (counts[1] + counts[3]))
      } else if (Reference[1] == "T") {
        return(counts[4] / (counts[4] + counts[2]))
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })]
}

# Add Uploaded_variation column
reditools[, Uploaded_variation := paste0(Region, "_", Position, 
                                         ifelse(Reference == "A", "_A/G", 
                                         ifelse(Reference == "T", "_T/C", "")))]

# Remove unnecessary columns
reditools <- reditools[, -c(1:324), with = FALSE]


##### RNA edited sites filtered with DNA mutations A/G and T/C ############################

# Load required libraries
library(dplyr)
library(tidyr)

# Define the working directory and folder path for MAF files
setwd("path_to_working_directory/")
folder_path <- "path_to_mutect2_vep_mafs/"
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


# Clinical data processing
deg <- read.delim("path_to_deg_file/v1_all_data_pivoted_filtered_104_subject_id.txt")
rownames(deg) <- deg$Name
deg$Name <- NULL
Run1 <- colnames(deg)

clinical_data_forR <- read.delim("path_to_clinical_data/clinical_data_forR.tsv")
BEAUTY_SraRunTable <- read.csv("path_to_BEAUTY_SraRunTable.csv")
phs001050.v1.p1 <- read.delim("path_to_phs001050_file.txt")

# Merge clinical data
merge <- merge(BEAUTY_SraRunTable, clinical_data_forR, by.x = "submitted_subject_id", by.y = "Subject_ID", all = TRUE)
merge_filtered <- subset(merge, submitted_subject_id %in% Run1) %>%
  filter(analyte_type == "DNA" & body_site == "Breast")
merge_filtered_4_reditools <- filter(merge, analyte_type == "RNA")

# Filter MAF for relevant samples
combined_maf <- subset(combined_maf, Sample %in% merge_filtered$Run)
maf_pass <- filter(combined_maf, FILTER == "PASS" & Variant_Type == "SNP")

# Apply specific filtering conditions
filtered_maf <- maf_pass %>%
  filter(
    (Reference_Allele == "A" & (Tumor_Seq_Allele1 == "G" | Tumor_Seq_Allele2 == "G")) |
      (Reference_Allele == "T" & (Tumor_Seq_Allele1 == "C" | Tumor_Seq_Allele2 == "C"))
  ) %>%
  mutate(
    Variant_ID = case_when(
      Reference_Allele != Tumor_Seq_Allele1 ~ paste0(Chromosome, "_", Start_Position, "_", Reference_Allele, "/", Tumor_Seq_Allele1),
      Reference_Allele != Tumor_Seq_Allele2 ~ paste0(Chromosome, "_", Start_Position, "_", Reference_Allele, "/", Tumor_Seq_Allele2),
      TRUE ~ paste0(Chromosome, "_", Start_Position, "_", Reference_Allele, "/", Tumor_Seq_Allele1)
    )
  )

write.table(filtered_maf, "filtered_maf_DNA_tumoral_editions.tsv", row.names = FALSE)

# Process germline MAF files
folder_path_germline <- "path_to_normal_sample_vep_mafs/"
maf_files_germline <- list.files(path = folder_path_germline, pattern = "\\.maf$", full.names = TRUE)

# Read and combine germline MAF files
maf_list_germline <- lapply(maf_files_germline, read_maf_file)
combined_maf_g <- do.call(rbind, maf_list_germline)
combined_maf_g$Sample <- gsub("\\.maf$", "", combined_maf_g$Sample)

# Filter germline MAF for specific criteria
maf_pass_g <- filter(combined_maf_g, FILTER == "PASS" & Variant_Type == "SNP")

filtered_maf_g <- maf_pass_g %>%
  filter(
    (Reference_Allele == "A" & (Match_Norm_Seq_Allele1 == "G" | Match_Norm_Seq_Allele2 == "G")) |
      (Reference_Allele == "T" & (Match_Norm_Seq_Allele1 == "C" | Match_Norm_Seq_Allele2 == "C"))
  ) %>%
  mutate(
    Variant_ID = case_when(
      Reference_Allele != Tumor_Seq_Allele1 ~ paste0(Chromosome, "_", Start_Position, "_", Reference_Allele, "/", Tumor_Seq_Allele1),
      Reference_Allele != Tumor_Seq_Allele2 ~ paste0(Chromosome, "_", Start_Position, "_", Reference_Allele, "/", Tumor_Seq_Allele2),
      TRUE ~ paste0(Chromosome, "_", Start_Position, "_", Reference_Allele, "/", Tumor_Seq_Allele1)
    )
  )

# Process reditools data
reditools <- read.delim("path_to_reditools_input/BEAUTY107_master_table_biallelic_minimal.tsv")
reditools <- filter(reditools, Reference == "A" | Reference == "T")

# Add BaseCount columns
base_count_columns <- grep("BaseCount.A.C.G.T.", names(reditools), value = TRUE)
for (sample_col in base_count_columns) {
  sample_name <- sub("BaseCount.A.C.G.T._", "", sample_col)
  reditools[[paste0("BaseCount_ref_", sample_name)]] <- sapply(seq_len(nrow(reditools)), function(i) {
    counts <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", reditools[i, sample_col]), ",")))
    if (length(counts) == 4) {
      ifelse(reditools$Reference[i] == "A", counts[1], ifelse(reditools$Reference[i] == "T", counts[4], NA))
    } else {
      NA
    }
  })
  reditools[[paste0("BaseCount_alt_", sample_name)]] <- sapply(seq_len(nrow(reditools)), function(i) {
    counts <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", reditools[i, sample_col]), ",")))
    if (length(counts) == 4) {
      ifelse(reditools$Reference[i] == "A", counts[3], ifelse(reditools$Reference[i] == "T", counts[2], NA))
    } else {
      NA
    }
  })
}

reditools <- reditools %>%
  mutate(Uploaded_variation = paste0(Region, "_", Position,
                                     ifelse(Reference == "A", "_A/G", ifelse(Reference == "T", "_T/C", "")))) %>%
  filter(!(Uploaded_variation %in% filtered_maf$Variant_ID))

# Rename and adjust SRR columns
srr_data <- reditools %>%
  select(all_of(c("Uploaded_variation", grep("SRR", names(reditools), value = TRUE))))

colnames(srr_data) <- gsub("BaseCount_(ref|alt)_", "", colnames(srr_data))
write.table(srr_data, "srr_data_ref_alt_reditools_filtered.tsv", row.names = FALSE)

# Calculate proportions and finalize dataset
for (srr_id in unique(gsub("BaseCount_(ref|alt)_", "", colnames(srr_data)[grepl("BaseCount_", colnames(srr_data))]))) {
  ref_col <- paste0("BaseCount_ref_", srr_id)
  alt_col <- paste0("BaseCount_alt_", srr_id)
  srr_data[[srr_id]] <- srr_data[[alt_col]] / (srr_data[[alt_col]] + srr_data[[ref_col]])
}




