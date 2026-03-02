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

cat("clinical_clean data with dimensions:", dim(clinical_clean), "\n")#1111 15
head(clinical_clean)

##--------------------------------------------------
## 3. Remove duplicate rows based on patient and barcode
##--------------------------------------------------

cat("Unique patients in clinical data:", length(unique(clinical_clean$patient)), "\n")
cat("Unique barcodes in clinical data:", length(unique(clinical_clean$barcode)), "\n")

table(table(clinical_clean$patient))

clinical_unique <- clinical_clean %>%
  distinct(patient, .keep_all = TRUE)
  
cat("clinical data of unique patients:", length(unique(clinical_unique$patient)), "\n")#1095 
cat("clinical data of unique barcodes:", length(unique(clinical_unique$barcode)), "\n")#1095 


##---------------------------------------------------
## 4. Missing data summary
##---------------------------------------------------

cat("Summary statistics of clinical_unique data:\n")
summary(clinical_unique)

cat("structure of clinical_unique data:\n")
str(clinical_unique)

#check the missing data in clinical_unique
cat("Missing data summary for clinical_unique:\n")

summarize_missing <- function(df){
  result <- df %>%
  summarise_all(~ sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_count") %>%
  mutate(missing_percentage = (missing_count / nrow(clinical_unique)) * 100) %>% 
  arrange(desc(missing_count))
  return(result)
}

missing_summary <-summarize_missing(clinical_unique)
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
cat("clinical data after filtering for missing values has dimensions:", dim(clinical_filtered_1), "\n")
cat ("Number of removed samples:", nrow(clinical_unique) - nrow(clinical_filtered), "\n")


missing_summary <-summarize_missing(clinical_filtered_1)
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
summary(clinical_filtered_2$days_to_last_follow_up)


#create prorper survival time and event variables
cat("\n vital status table for clinical_filtered_2:\n")
table(clinical_filtered_2$vital_status)

clinical_filtered_2 <- clinical_filtered_2 %>%
  mutate(
    # Overall survival time (in days)
    OS_time = case_when(
       vital_status == "Dead" & !is.na(days_to_death) ~ days_to_death,
       vital_status == "Alive" & !is.na(days_to_last_follow_up)  ~ days_to_last_follow_up,
      TRUE ~ NA_real_
    ),
    # Overall survival event (1=death, 0=censored/alive)
    OS_event = case_when(
      vital_status == "Alive" ~ 0,
      vital_status == "Dead" ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Convert to years and months for interpretability
    OS_time_years = OS_time / 365.25,
    OS_time_months = OS_time / 30.44
  )

cat("dimension of clinical_filtered_2:", dim(clinical_filtered_2), "\n")#1036 19 

head(clinical_filtered_2)
summary(clinical_filtered_2$OS_time)
table(clinical_filtered_2$OS_time, useNA="ifany")

# Remove samples with invalid/missing survival data and remove negative or zero survival times
clinical_filtered_3 <- clinical_filtered_2 %>%
  filter(!is.na(OS_time) & OS_time > 0 )

cat("dimension of clinical_filtered_3:", dim(clinical_filtered_3), "\n")#1014 19  



cat("\nmissing_summary of clinical_filtered_3:\n")

missing_summary <-summarize_missing(clinical_filtered_3)
print(missing_summary)


##---------------------------------------------------------------------
## 7. Clean and standardize PAM50 Subtype variable
##---------------------------------------------------------------------

cat("\nCleaning PAM50 Subtype variable:\n")

# Check unique PAM50 values
cat("Unique PAM50_Subtype values before cleaning:\n")
unique(clinical_filtered_3$PAM50_Subtype)

clinical_filtered_3 %>%
  count(PAM50_Subtype, .drop=FALSE)



##---------------------------------------------------------------------
## 8. Clean and standardize paper pathologic_stage variable
##---------------------------------------------------------------------
cat("\nCleaning pathologic_stage variable:\n")

# # Check unique pathologic_stage values
# cat("Unique pathologic_stage values before cleaning:\n")
# unique(clinical_filtered_3$pathologic_stage)

# clinical_filtered_3 %>%
#   count(pathologic_stage, .drop=FALSE)
  
cat("Unique paper pathologic_stage values before cleaning:\n")
unique(clinical_filtered_3$paper_pathologic_stage)

clinical_filtered_3 %>%
  count(paper_pathologic_stage, .drop=FALSE)

#replace NA with unknown in paper_pathologic_stage

replace_NA_unknown <- function(df, column_name) {
  
  # Convert string column name to a symbol for tidy evaluation
  col_sym <- sym(column_name)
  
  result <- df %>%
    mutate(!!col_sym := case_when(
      is.na(!!col_sym) | !!col_sym == "NA" ~ "unknown",
      TRUE ~ !!col_sym
    ))
  
  # Print the count
  cat("Counts for", column_name, ":\n")
  print(result %>% count(!!col_sym, .drop = FALSE))
  
  return(result)
}
clinical_filtered_3 <- replace_NA_unknown(clinical_filtered_3, "paper_pathologic_stage")



##---------------------------------------------------------------------
## 9. Clean and standardize TNM variable
##---------------------------------------------------------------------

#ajcc_pathologic_t
cat("\nCleaning TNM variables:\n")

cat("Unique ajcc_pathologic_t values before cleaning:\n")
table(clinical_filtered_3$ajcc_pathologic_t, useNA="ifany")

clinical_filtered_3 <- replace_NA_unknown(clinical_filtered_3, "ajcc_pathologic_t")


#ajcc_pathologic_n

cat("Unique ajcc_pathologic_n values before cleaning:\n")
table(clinical_filtered_3$ajcc_pathologic_n, useNA="ifany")

clinical_filtered_3 <- replace_NA_unknown(clinical_filtered_3, "ajcc_pathologic_n")


#ajcc_pathologic_m
cat("Unique ajcc_pathologic_m values before cleaning:\n")
table(clinical_filtered_3$ajcc_pathologic_m, useNA="ifany")

clinical_filtered_3 <- replace_NA_unknown(clinical_filtered_3, "ajcc_pathologic_m")


##---------------------------------------------------------------------
## 10. Remove uncessary variables
##---------------------------------------------------------------------
clinical_filtered_4 <- clinical_filtered_3 %>%
  select(-c("pathologic_stage", "sample_type"))
  
cat("clinical_filtered_4 dimensions:", dim(clinical_filtered_4), "\n")#1034 17

saveRDS(clinical_filtered_4, "data/mid_files/clinical_filtered_4.rds")


##---------------------------------------------------------------------
## 11. Create age group
##---------------------------------------------------------------------
clinical_filtered_4 <-clinical_filtered_4 %>%
  mutate(
    age_at_diagnosis_years =age_at_diagnosis/365.25,
    age_group= case_when(
      age_at_diagnosis_years <40 ~ "<40",
      age_at_diagnosis_years <50 ~ "40-49",
      age_at_diagnosis_years <60 ~ "50-59",
      age_at_diagnosis_years <70 ~ "60-69",
      TRUE ~"70+" 
    )
  )
head(clinical_filtered_4)

table(clinical_filtered_4$age_group, useNA="always")

#final summary
summary(clinical_filtered_4)

##---------------------------------------------------------------------
## 12. Extract the matching expression data (unstarnded and tpm_unstrand)
##---------------------------------------------------------------------
#unstarnded
cat("Extract expression count data:\n")
counts_raw <- assay(BRCA_rna_data, "unstranded")
cat("Dimension of count_raw:", dim(counts_raw), "\n")

counts_match <- counts_raw [, clinical_filtered_4$barcode]
cat("Dimension of count_match:", dim(counts_match), "\n")

#tpm_unstrand
cat("tpm_unstrand:\n")
tpm_unstrand_raw <- assay(BRCA_rna_data,"tpm_unstrand")
cat("tpm_unstrand_raw:", dim(tpm_unstrand_raw), "\n")


tpm_unstrand_match <- tpm_unstrand_raw [, clinical_filtered_4$barcode]
cat("Dimension of tpm_unstrand_match:", dim(tpm_unstrand_match), "\n")


##---------------------------------------------------------------------
## Save final filtered clinical data and expression matrix
##---------------------------------------------------------------------
saveRDS(clinical_filtered_4, "data/processed/clinical_filtered_clean.rds")
saveRDS(counts_match, "data/processed/counts_match.rds")
saveRDS(tpm_unstrand_match, "data/processed/tpm_unstrand_match.rds")





