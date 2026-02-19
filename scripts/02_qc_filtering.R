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



## --------------------------------------------------------
## 2. Extract, process and quality control for clinical data 
## ----------------------------------------------------------
clinical <- colData(BRCA_rna_data) %>%
  as.data.frame() 
cat("Extracted clinical data with dimensions:", dim(clinical), "\n")
colnames(clinical)


list_cols <- sapply(clinical, is.list)
clinical_1 <- clinical[, !list_cols]

write.csv(clinical_1, "data/processed/clinical_1.csv", row.names = FALSE)

##------------------------------------------------
## 5a. Extract key variables
##------------------------------------------------
clinical_clean <- clinical %>%
  select(
    patient = patient,
    barcode = barcode,
    sample_type = sample_type,
    age_at_diagnosis = age_at_diagnosis,
    race= race,
    gender=gender,
    country_residence = country_of_residence_at_enrollment,
    vital_status = vital_status,
    days_to_death = days_to_death,
    days_to_last_follow_up = days_to_last_follow_up, #patient alive
    paper_pathologic_stage=paper_pathologic_stage,
    pathologic_stage = ajcc_pathologic_stage,
    ajcc_pathologic_t = ajcc_pathologic_t,
    ajcc_pathologic_n = ajcc_pathologic_n,
    ajcc_pathologic_m = ajcc_pathologic_m,
    PAM50_Subtype = paper_BRCA_Subtype_PAM50
  )

cat("clinical_clean data with dimensions:", dim(clinical_clean), "\n")
head(clinical_clean)

##--------------------------------------------------
## 5b. Remove duplicate rows based on patient and barcode
##--------------------------------------------------

cat("Unique patients in clinical data:", length(unique(clinical_clean$patient)), "\n")
cat("Unique barcodes in clinical data:", length(unique(clinical_clean$barcode)), "\n")

clinical_clean <- clinical_clean %>%
  distinct(patient, barcode, .keep_all = TRUE)




##---------------------------------------------------
##check the distribution and levels of key variables 
##---------------------------------------------------

#pathologic_stage
cat("pathologic_stage distribution:\n")
clinical_clean %>%
           count(pathologic_stage, .drop=FALSE) 
 
str(clinical_clean$pathologic_stage)

#PAM50_Subtype
cat("PAM50_Subtype distribution:\n")
clinical_clean %>%
           count(PAM50_Subtype, .drop=FALSE) 

#primary_diagnosis
cat("primary diagnosis distribution:\n")
clinical_clean %>%
           count(primary_diagnosis, .drop=FALSE)





table(clinical$tumor_descriptor)

str(clinical$tumor_descriptor)


table(clinical_clean$paper_BRCA_Subtype_PAM50)
table(clinical_clean$ajcc_pathologic_stage)
table(clinical_clean$vital_status)







  select(-c("barcode", "sample_type", "project_id", "data_category", "data_type", "platform")) %>%
  distinct() # Remove duplicate rows