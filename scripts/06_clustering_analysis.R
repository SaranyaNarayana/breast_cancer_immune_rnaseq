################################################################################
# 04_clustering_analysis.R
# Goal: Identify immune subtypes within each PAM50 breast cancer subtype
#
# Strategy:
#   - Cluster SAMPLES by their immune feature profiles (not the features)
#   - Use PCA reduction first to remove correlation between immune features
#   - Run ConsensusClusterPlus on PC scores for stability
#   - Where clustering fails (imbalanced), use continuous immune scores
#
# Output:
#   - Discrete immune clusters (IC1/IC2) for Her2 and Normal subtypes
#   - Continuous immune scores + binary groups for Basal, LumA, LumB
#   - Unified group column for downstream survival + characterization
################################################################################


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
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(cowplot)
  library(ggpubr)
  library(tibble)
  library(scales)
  library(ggrepel)
})


# Create output directories
dir.create("results/figures/clustering",       recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/clustering/pca",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/clustering/umap",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/consensus_clustering",     recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/clustering",        recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/clustering",        recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# 1. LOAD DATA
# Why: Need immune features merged with clinical metadata (PAM50 subtype)
#      CIBERSORT has already been separated out in script 03
#===============================================================================

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

head(data_full[, 155:170])

cat("dimension of data full:", dim(data_full), "\n")#1013 178

cat("PAM50 Subtype distribution:\n")
table(data_full$PAM50_Subtype)

saveRDS(data_full, "data/processed/clustering/data_full.rds")


## ---------------------------------------------------------------------------
## 2. CIBERSORT  QC filtering and used for validation of immune-based clusters
## ---------------------------------------------------------------------------


cibersort_col <- grep("^CIBERSORT_", colnames(data_full), value=TRUE)

cibersort_immune_estimates <- data_full %>%
    select(sample, all_of(cibersort_col))
cat("Dimension of cibersort_immune_estimates:", dim(cibersort_immune_estimates), "\n")


# Apply QC filter and extract validation dataset for downstream analysis
cibersort_immune_estimates %>%
    select(`CIBERSORT_P-value`, `CIBERSORT_RMSE`,  `CIBERSORT_Correlation`) %>%
    summary()


cibersort_valid <- cibersort_immune_estimates %>%
  filter(`CIBERSORT_P-value` < 0.05,
         `CIBERSORT_Correlation` > 0,
         `CIBERSORT_RMSE` < quantile(`CIBERSORT_RMSE`, 0.95, na.rm = TRUE))
cat("After CIBERSORT QC filtering, retained", dim(cibersort_valid),"\n") # 670 26 

saveRDS(cibersort_valid, "data/processed/clustering/cibersort_valid.rds")


## ---------------------------------------------------------------------------
# 3. Remove ALL CIBERSORT columns from main data
## ---------------------------------------------------------------------------
data_main <- data_full %>%
  select(-starts_with("CIBERSORT_"))
cat("Dimension of data main after removing CIBERSORT columns:", dim(data_main), "\n")#1013 153

        

## ---------------------------------------------------------------------------
## 4. Define clustering columns (immune features only, no clinical metadata)
## ---------------------------------------------------------------------------

clustering_feature_cols <- grep(
  "^(GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_)",
  colnames(data_main), value = TRUE
)
cat("Number of immune features used for clustering:", length(clustering_feature_cols), "\n")# 134
