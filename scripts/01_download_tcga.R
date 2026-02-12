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
## 1. Install and load libraries
## -----------------------------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
install.packages(c("dplyr", "tibble"))

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
# cat("===== GDC download =====\n")
# log_message("Starting GDC download")

# GDCdownload(
#   query,
#   method = "api",
#   files.per.chunk = 20
# )

# log_message("Finished GDC download")


#if downloaded already use the below code to prepare the data
GDCdownload(query, method = "api", files.per.chunk = 20, directory = "GDCdata/")

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

  saveRDS(
    se_batch,
    file = paste0("data/raw/brca_se_batch_", i, ".rds")
  )

  #CRITICAL PART_Free meemory after each batch to avoid OOM errors
  rm(se_batch, query_batch)
  gc()

  cat("Finished batch", i, "\n")
}


#merging progressively due to memory constraints, we will read in the batch files and combine them into one SE object
files <- list.files("data/raw", pattern = "brca_se_batch_.*rds", full.names = TRUE)
length(files) # check number of batch files

files <- sort(files)

meta_names <- lapply(files, function(f) {
  se <- readRDS(f)
  cn <- colnames(colData(se))
  rm(se)
  gc()
  return(cn)
})


# Get list of metadata data.frames:there difference in the metadata columns across batches, so we will need to handle that when combining SE objects
meta_list <- lapply(se_list, function(se) as.data.frame(colData(se)))

# Get intersection of column names
common_cols <- Reduce(intersect, lapply(meta_list, colnames))
length(common_cols)#85

#Build Metadata Incrementally
meta_full <- NULL

for (f in files) {
  cat("Processing:", f, "\n")
  
  se <- readRDS(f)
  meta <- as.data.frame(colData(se))
  
  # Keep only common columns safely
  meta <- meta[, intersect(colnames(meta), common_cols), drop = FALSE]
  
  meta_full <- rbind(meta_full, meta)
  
  rm(se, meta)
  gc()
}

dim(meta_full)#1111   76
saveRDS(meta_full, "data/raw/meta_full.rds")



#Combine Expression Matrix
# Get common genes first
gene_lists <- lapply(files, function(f) {
  se <- readRDS(f)
  g <- rownames(se)
  rm(se)
  gc()
  return(g)
})

common_genes <- Reduce(intersect, gene_lists)
length(common_genes)#60660

expr_full <- NULL

for (f in files) {
  cat("Processing:", f, "\n")
  
  se <- readRDS(f)
  mat <- assay(se)
  
  mat <- mat[common_genes, , drop = FALSE]
  
  expr_full <- cbind(expr_full, mat)
  
  rm(se, mat)
  gc()
}

dim(expr_full)#60660  1111
saveRDS(expr_full, "data/raw/expr_full.rds")


#Keep genes with counts ≥10 in at least 10% of samples

keep_genes <- rowSums(expr_full >= 10) >= (0.10 * ncol(expr_full))

expr_filtered <- expr_full[keep_genes, ]

dim(expr_filtered)#26006  1111

#Rebuild Clean SummarizedExperiment

brca_se_clean <- SummarizedExperiment(
  assays = list(counts = expr_filtered),
  colData = meta_full
)

brca_se_clean
saveRDS(brca_se_clean, "data/processed/brca_se_clean.rds")




















































































# brca_se <- GDCprepare(query)

# ## Save raw object for reproducibility
# saveRDS(brca_se, "data/raw/TCGA_BRCA_STAR_counts_SE.rds")

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
