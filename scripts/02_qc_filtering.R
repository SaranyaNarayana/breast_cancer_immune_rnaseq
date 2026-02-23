############################################################
# Quality filtering and preprocessing of TCGA-BRCA RNA-seq data
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Derive key variables from clinical data
#   - Remove duplicate samples and ensure unique patient representation
#   - Check distribution of key clinical variables and remove missing data if necessary



#   - Quality control and filtering of TCGA-BRCA RNA-seq data
#   - Prepare clean count matrix and metadata for downstream analysis
#   - Save processed data for analysis
#   - Conduct exploratory data analysis (EDA) and visualization
############################################################

############################################################
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
## 3. Extract key clinical variables
##------------------------------------------------
clinical_clean <- clinical %>%
  select(
    patient = patient,
    barcode = barcode,
    sample_type = sample_type,
    age_at_diagnosis = age_at_diagnosis,
    race= race,
    gender=gender,
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
## 3. Remove duplicate rows based on patient and barcode
##--------------------------------------------------

cat("Unique patients in clinical data:", length(unique(clinical_clean$patient)), "\n")
cat("Unique barcodes in clinical data:", length(unique(clinical_clean$barcode)), "\n")

table(table(clinical_clean$patient))

clinical_unique <- clinical_clean %>%
  distinct(patient, .keep_all = TRUE)
  
cat("clinical data of unique patients:", length(unique(clinical_unique$patient)), "\n")
cat("clinical data of unique barcodes:", length(unique(clinical_unique$barcode)), "\n")


##---------------------------------------------------
## 4. Missing data summary
##---------------------------------------------------

cat("Summary statistics of clinical_unique data:\n")
summary(clinical_unique)

cat("structure of clinical_unique data:\n")
str(clinical_unique)

#check the missing data in clinical_unique
cat("Missing data summary for clinical_unique:\n")
missing_summary <- clinical_unique %>%
  summarise_all(~ sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_count") %>%
  mutate(missing_percentage = (missing_count / nrow(clinical_unique)) * 100) %>% 
  arrange(desc(missing_count))

print(missing_summary)




##---------------------------------------------------------------------
## 5. Remove smaples with missing critical variables 
##(vital status, gender, age_at_diagnosis and PAM50 subtype; All 4 vairable must be present for a sample to be retained)
##---------------------------------------------------------------------
cat("\nHandling missing values:\n")

clinical_filtered_1 <- clinical_unique %>%
  filter(
    !is.na(vital_status) & vital_status != "",
    !is.na(gender) & gender != "",
    !is.na(age_at_diagnosis),
    !is.na(PAM50_Subtype) & PAM50_Subtype != ""
  )
cat("clinical data after filtering for missing values has dimensions:", dim(clinical_filtered), "\n")
cat ("Number of removed samples:", nrow(clinical_unique) - nrow(clinical_filtered), "\n")


missing_summary <- clinical_filtered_1  %>%
  summarise_all(~ sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_count") %>%
  mutate(missing_percentage = (missing_count / nrow(clinical_unique)) * 100) %>% 
  arrange(desc(missing_count))

print(missing_summary)





##---------------------------------------------------------------------
## 6. Handle survival variables
##---------------------------------------------------------------------
cat("\nHandling survival variables:\n")

#Fix mssing values in days_to_last_follow_up using days_to_death for patients who are alive but have missing follow-up time
clinical_filtered_2 <- clinical_filtered_1 %>%
  mutate(
    days_to_last_follow_up=case_when(
      is.na(days_to_last_follow_up) &
      vital_status == "Alive" &
      !is.na(days_to_death) ~ days_to_death,
      TRUE ~ days_to_last_follow_up   
    )
  )
saveRDS(clinical_filtered_2, "data/mid_files/clinical_filtered_2.rds")


#create prorper survivaal time and event variables























saveRDS(unstranded_full, "data/mid_files/unstranded_full.rds")





















#pathologic_stage
cat("pathologic_stage distribution:\n")
clinical_unique %>%
           count(pathologic_stage, .drop=FALSE) 
 
str(clinical_unique$pathologic_stage)

#PAM50_Subtype
cat("PAM50_Subtype distribution:\n")
clinical_unique %>%
           count(PAM50_Subtype, .drop=FALSE) 






table(clinical$tumor_descriptor)

str(clinical$tumor_descriptor)


table(clinical_clean$paper_BRCA_Subtype_PAM50)
table(clinical_clean$ajcc_pathologic_stage)
table(clinical_clean$vital_status)


##--------------------------------------------------
## 5d. PAM50 Subtype distribution across pathologic stages
##--------------------------------------------------




  select(-c("barcode", "sample_type", "project_id", "data_category", "data_type", "platform")) %>%
  distinct() # Remove duplicate rows


#Then subset expression matrix:
  expr_unique <- expr[, clinical_unique$barcode]