############################################################
# Identify immune-based subtypes within each clinical subtype 
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Define Immune signatures from MSigDB and literature
#   - Convert Ensembl IDs to gene symbols
#   - Run gene set variation analysis (GSVA) and single sample GSEA (ssGSEA) for immune signatures
#   - Run immune deconvolution using quanTIseq, MCP-counter, EPIC, xCell, and CIBERSORT 
#   - Combine all resutls for further analysis and visualization in downstream scripts



############################################################
setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")
############################################################

#-----------------------------------
# Logging
#-----------------------------------

log_file <- file(
  paste0("logs/06_cluster_analysis_",
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
  library(ConsensusClusterPlus)
  library(cluster)
  library(factoextra)
  library(umap)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(cowplot)
})

install.packages("cowplot")



#create output directory
dir.create("results/figures/clustering", recursive=TRUE, showWarnings = FALSE)
dir.create("data/processed/clustering", recursive=TRUE, showWarnings = FALSE)
dir.create("results/consensus_clustering", recursive=TRUE, showWarnings = FALSE)


## -----------------------------
## 1. Load Data
## -----------------------------

immune_features <- readRDS("data/processed/immune/immune_features_combined.rds")
cat("Loaded immune features data with dimensions:", dim(immune_features), "\n")#1013 160
head(immune_features[, 155:160])
colnames(immune_features)


clinical <- readRDS("data/processed/clinical_filtered.rds")
cat("Loaded clinical data with dimensions:", dim(clinical), "\n")#1013 19
head(clinical[, 1:5])
colnames(clinical)

# merge immune features with clinical data
data_full <- immune_features %>%
  left_join(clinical, by = c("sample" = "barcode"))

cat("dimension of data full:", dim(data_full), "\n")

table(data_full$PAM50_Subtype)
