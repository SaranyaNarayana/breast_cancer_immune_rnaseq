############################################################
# Quality filtering, preprocessing and plotting of clinical data 
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Derive key variables from clinical data and conduct exploratory data analysis (EDA) and visualization
#   - Remove duplicate samples and ensure unique patient representation
#   - Remove missing data if necessary



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
  library(ggplot2)
  library(patchwork)  # For combining plots
  library(RColorBrewer)
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

#function to replace NA with unknown in paper_pathologic_stage

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
  
cat("clinical_filtered_4 dimensions:", dim(clinical_filtered_4), "\n")#1014 17 

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
## 13. Save final filtered clinical data and expression matrix
##---------------------------------------------------------------------
saveRDS(clinical_filtered_4, "data/processed/clinical_filtered_clean.rds")
saveRDS(counts_match, "data/processed/counts_match.rds")
saveRDS(tpm_unstrand_match, "data/processed/tpm_unstrand_match.rds")



##---------------------------------------------------------------------
## 14. Final plots for clinical data
##---------------------------------------------------------------------
colnames(clinical_filtered_4)
dir.create("results/figures/clinical_overview", recursive = TRUE, showWarnings = FALSE)

##---------------------------------------------------------
## Plot 1: DEMOGRAPHIC OVERVIEW-Age distribution
##---------------------------------------------------------

p1 <- ggplot(clinical_filtered_4, aes(x = age_at_diagnosis_years)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(age_at_diagnosis_years)), 
             color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = median(clinical_filtered_4$age_at_diagnosis_years) + 5, 
           y = Inf, vjust = 2,
           label = paste("Median:", round(median(clinical_filtered_4$age_at_diagnosis_years), 1), "years"),
           color = "red", size = 4) +
  labs(title = "Age Distribution at Diagnosis",
       x = "Age (years)", y = "Number of Patients") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/01_age_distribution.png", 
       p1, width = 8, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 2: DEMOGRAPHIC OVERVIEW- Age by PAM50 subtype
##---------------------------------------------------------

p2 <- ggplot(clinical_filtered_4, aes(x = PAM50_Subtype, y = age_at_diagnosis_years, 
                                        fill = PAM50_Subtype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1)),
               vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Age Distribution by PAM50 Subtype",
       x = "PAM50 Subtype", y = "Age at Diagnosis (years)") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/02_age_by_pam50.png", 
       p2, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 3: DEMOGRAPHIC OVERVIEW- Age groups distribution
##---------------------------------------------------------
p3 <- ggplot(clinical_filtered_4, aes(x = age_group, fill = age_group)) +
  geom_bar(alpha = 0.8) +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "Distribution by Age Group",
       x = "Age Group", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/03_age_groups.png", 
       p3, width = 8, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 4: PAM50 SUBTYPE ANALYSIS:PAM50 subtype distribution
##---------------------------------------------------------
pam50_counts <- clinical_filtered_4 %>%
  count(PAM50_Subtype) %>%
  mutate(percentage = n / sum(n) * 100)

p4 <- ggplot(pam50_counts, aes(x = reorder(PAM50_Subtype, -n), y = n, 
                                 fill = PAM50_Subtype)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "PAM50 Molecular Subtype Distribution",
       subtitle = paste0("Total patients: ", nrow(clinical_filtered_4)),
       x = "PAM50 Subtype", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/04_pam50_distribution.png", 
       p4, width = 10, height = 6, dpi = 300)

##---------------------------------------------------------
## Plot 5: PAM50 SUBTYPE ANALYSIS:PAM50 pie chart
##---------------------------------------------------------  

p5 <- ggplot(pam50_counts, aes(x = "", y = n, fill = PAM50_Subtype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(PAM50_Subtype, "\n", round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "PAM50 Subtype Distribution (Pie Chart)") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "right")

ggsave("results/figures/clinical_overview/05_pam50_pie.png", 
       p5, width = 10, height = 8, dpi = 300)


##---------------------------------------------------------
## Plot 6: STAGE ANALYSIS:Stage distribution
##---------------------------------------------------------

stage_counts <- clinical_filtered_4 %>%
  count(paper_pathologic_stage) %>%
  mutate(percentage = n / sum(n) * 100)

p6 <- ggplot(stage_counts, aes(x = reorder(paper_pathologic_stage, -n), 
                                 y = n, fill = paper_pathologic_stage)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Pathologic Stage Distribution",
       x = "Pathologic Stage", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/clinical_overview/06_stage_distribution.png", 
       p6, width = 10, height = 6, dpi = 300)

##---------------------------------------------------------
## Plot 7: STAGE ANALYSIS:Stage by PAM50 subtype
##---------------------------------------------------------
p7 <- ggplot(clinical_filtered_4, aes(x = PAM50_Subtype, fill = paper_pathologic_stage)) +
  geom_bar(position = "fill", alpha = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Stage Distribution by PAM50 Subtype",
       x = "PAM50 Subtype", y = "Proportion of Patients",
       fill = "Stage") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/07_stage_by_pam50.png", 
       p7, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 8: STAGE ANALYSIS:Stage by PAM50 subtype
##---------------------------------------------------------

p8 <- ggplot(clinical_filtered_4, aes(x = PAM50_Subtype, fill = paper_pathologic_stage)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Stage Distribution by PAM50 Subtype (Counts)",
       x = "PAM50 Subtype", y = "Number of Patients",
       fill = "Stage") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/08_stage_by_pam50_counts.png", 
       p8, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 9: Violin plot: PAM50 subtype vs age at diagnosis
##---------------------------------------------------------

p9 <-ggplot(clinical_filtered_4, 
       aes(x = PAM50_Subtype, y = age_at_diagnosis_years, fill = PAM50_Subtype)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
  geom_jitter(alpha = 0.2, width = 0.1, size = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Age Distribution by PAM50 Subtype",
       x = "PAM50 Subtype", y = "Age at Diagnosis (years)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results/figures/clinical_overview/09_age_violin_pam50.png", p9, width = 10, height = 6, dpi = 300)

##---------------------------------------------------------
## Plot 10: Violin plot: stage vs age at diagnosis
##---------------------------------------------------------

p10 <- ggplot(clinical_filtered_4 %>% filter(paper_pathologic_stage != "Unknown"), 
       aes(x = paper_pathologic_stage, y = age_at_diagnosis_years, 
           fill = paper_pathologic_stage)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
  geom_jitter(alpha = 0.2, width = 0.1, size = 0.5) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Age Distribution by Stage",
       x = "Stage", y = "Age at Diagnosis (years)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results/figures/clinical_overview/10_age_violin_stage.png", p10, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 11: DEMOGRAPHIC CHARACTERISTICS - Gender distribution
##---------------------------------------------------------

p11 <- ggplot(clinical_filtered_4, aes(x = gender, fill = gender)) +
  geom_bar(alpha = 0.8) +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Female" = "#FF69B4", "Male" = "#4169E1")) +
  labs(title = "Gender Distribution",
       x = "Gender", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/11_gender_distribution.png", 
       p11, width = 8, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 11: DEMOGRAPHIC CHARACTERISTICS - Race distribution
##---------------------------------------------------------

race_counts <- clinical_filtered_4 %>%
  count(race) %>%
  arrange(desc(n))

p12 <- ggplot(race_counts, aes(x = reorder(race, n), y = n, fill = race)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = n), hjust = -0.3, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Race/Ethnicity Distribution",
       x = "Race", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/12_race_distribution.png", 
       p12, width = 10, height = 6, dpi = 300)
