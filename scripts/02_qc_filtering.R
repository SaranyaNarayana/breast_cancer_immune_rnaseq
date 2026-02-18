############################################################
# Quality filtering and preprocessing of TCGA-BRCA RNA-seq data
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Quality control and filtering of TCGA-BRCA RNA-seq data
#   - Prepare clean count matrix and metadata for downstream analysis
#   - Save processed data for analysis
#   - Conduct exploratory data analysis (EDA) and visualization
############################################################

############################################################
getwd()
setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")

############################################################

#-----------------------------------
# Logging
#-----------------------------------

log_file <- file(
  paste0("logs/02_qc_filtering_",
         format(Sys.time(), "%Y%m%d_%H%M%S"),
         ".log"),
  open = "wt"
)

sink(log_file, type = "output")
sink(log_file, type = "message")

log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}


suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
})




## -----------------------------
## 1. Load Data
## -----------------------------
BRCA_rna_data <- readRDS("data/processed/TCGA_BRCA_rna_data.rds")
cat("Loaded clean SummarizedExperiment with dimensions:", dim(BRCA_rna_data), "\n")


# meta_full <- readRDS("data/raw/meta_full.rds")
# cat("Loaded metadata with dimensions:", dim(meta_full), "\n")



## -----------------------------
## 2. Extract and process clinical data 
## -----------------------------
clinical <- colData(BRCA_rna_data) %>%
  as.data.frame() 
cat("Extracted clinical data with dimensions:", dim(clinical), "\n")
colnames(clinical)

# Extract key variables
clinical_clean <- clinical %>%
  dplyr::select(
    patient = patient,
    barcode = barcode,
    sample_type = sample_type,
    age_at_diagnosis = age_at_diagnosis,
    vital_status = vital_status,
    days_to_death = days_to_death,
    days_to_last_follow_up = days_to_last_follow_up,
    tumor_stage = tumor_stage,
    pathologic_stage = ajcc_pathologic_stage,
    grade = neoplasm_histologic_grade,
    ER_status = paper_BRCA_Subtype_PAM50,  # Adjust based on actual column names
    PR_status = paper_PR_Status,
    HER2_status = paper_HER2_Final_Status,
    PAM50 = paper_BRCA_Subtype_PAM50
  )

cat("Processed clinical data with dimensions:", dim(clinical_clean), "\n")




  select(-c("barcode", "sample_type", "project_id", "data_category", "data_type", "platform")) %>%
  distinct() # Remove duplicate rows