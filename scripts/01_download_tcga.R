############################################################
# TCGA BRCA RNA-seq download & preprocessing
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Download TCGA-BRCA RNA-seq data
#   - Prepare count matrix and metadata
#   - Filter low-expression genes
############################################################
getwd()
setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")

############################################################

#-----------------------------------
# Logging
#-----------------------------------
dir.create("logs", showWarnings = FALSE)

log_file <- file(
  paste0("logs/01_tcga_brca_download_",
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


## -----------------------------
## 1. Load libraries
## -----------------------------
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
})

## -----------------------------
## 2. Define output directories
## -----------------------------
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 3. Query TCGA BRCA RNA-seq
## -----------------------------
cat("===== gdc project information =====\n")
gdcproj <- getGDCprojects()  #check the available projects
cat("Data dimensions:")
dim(gdcproj)#91 10
head(gdcproj)

cat("===== TCGA-BRCA project information =====\n")
brca_summary <- getProjectSummary("TCGA-BRCA") # get project summary
brca_summary


query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)


## -----------------------------
## 4. Download data
## -----------------------------
cat("===== GDC download =====\n")
log_message("Starting GDC download")

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 20
)

log_message("Finished GDC download")


#if downloaded already use the below code to prepare the data
#GDCdownload(query, method = "api", files.per.chunk = 20, directory = "GDCdata/")

total_files <-nrow(query$results[[1]])
cat("Total number of files:", total_files, "\n")

cat("Unique patient samples:",length(unique(query$results[[1]]$cases)), "\n")


## -----------------------------
## 5. Prepare SummarizedExperiment
## -----------------------------

#Processing 200 samples batchwise to avoid memory issues

batch_size <- 200
n_batches <- ceiling(total_files / batch_size)

for (i in 1:n_batches) {
  
  cat("Processing batch", i, "of", n_batches, "\n")
  
  start_idx <- ((i - 1) * batch_size) + 1
  end_idx <- min(i * batch_size, total_files)
  
  query_batch <- query
  query_batch$results[[1]] <- query$results[[1]][start_idx:end_idx, ]
  
  se_batch <- GDCprepare(query_batch)
  
  saveRDS(se_batch, file = paste0("data/raw/brca_se_batch_", i, ".rds"))
}























brca_se <- GDCprepare(query)

## Save raw object for reproducibility
saveRDS(brca_se, "data/raw/TCGA_BRCA_STAR_counts_SE.rds")

## -----------------------------
## 6. Extract count matrix
## -----------------------------
counts <- assay(brca_se)

## Clean gene IDs (remove version numbers)
rownames(counts) <- gsub("\\..*", "", rownames(counts))

## -----------------------------
## 7. Extract metadata
## -----------------------------
metadata <- colData(brca_se) |> 
  as.data.frame() |> 
  rownames_to_column("sample_id")

## Keep clinically relevant variables
metadata_clean <- metadata |> 
  select(
    sample_id,
    patient = patient,
    gender,
    age_at_diagnosis,
    tumor_stage,
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    er_status_by_ihc,
    pr_status_by_ihc,
    her2_status_by_ihc
  )

## -----------------------------
## 8. Basic QC: remove low-count genes
## -----------------------------
keep_genes <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep_genes, ]

## -----------------------------
## 9. Save processed data
## -----------------------------
saveRDS(counts_filtered, "data/processed/TCGA_BRCA_counts_filtered.rds")
write.csv(metadata_clean, "data/processed/TCGA_BRCA_metadata.csv", row.names = FALSE)

## -----------------------------
## 10. Summary output
## -----------------------------
message("TCGA BRCA preprocessing complete")
message("Samples: ", ncol(counts_filtered))
message("Genes retained: ", nrow(counts_filtered))

## -----------------------------
## Close log
## -----------------------------
sink(type = "message")
sink(type = "output")
close(log_file)
