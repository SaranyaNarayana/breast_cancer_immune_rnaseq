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
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(cowplot)
})



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

head(data_full[, 155:170])

cat("dimension of data full:", dim(data_full), "\n")#1013 178

cat("PAM50 Subtype distribution:\n")
table(data_full$PAM50_Subtype)


## ---------------------------------------------------------------------------
## 2. CIBERSORT  QC filtering and used for validation of immune-based clusters
## ---------------------------------------------------------------------------
cat("Performing CIBERSORT QC filtering...\n")

# Visualize CIBERSORT QC metrics
p_pval <- ggplot(data_full, aes(x = `CIBERSORT_P-value`)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "CIBERSORT P-value Distribution",
       x = "P-value", y = "Count") +
  theme_bw()

p_cor <- ggplot(data_full, aes(x = `CIBERSORT_Correlation`)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "CIBERSORT Correlation Distribution",
       x = "Correlation", y = "Count") +
  theme_bw()

p_rmse <- ggplot(data_full, aes(x = `CIBERSORT_RMSE`)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  geom_vline(xintercept = quantile(data_full$`CIBERSORT_RMSE`, 0.95, na.rm = TRUE),
             color = "red", linetype = "dashed") +
  labs(title = "CIBERSORT RMSE Distribution",
       x = "RMSE", y = "Count") +
  theme_bw()

p_qc <- plot_grid(p_pval, p_cor, p_rmse, ncol = 3)
ggsave("results/figures/clustering/CIBERSORT_QC_metrics.png",
       p_qc, width = 15, height = 5)

  

# Extract CIBERSORT data seperately 

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
## 3. Feature selection for clustering analysis 
## Using: GSVA, ssGSEA, quanTIseq, MCP, EPIC, Xcell only for clustering analysis
## ---------------------------------------------------------------------------

immune_features_for_clustering <- data_full %>%
  select(sample, starts_with("GSVA_"), starts_with("ssGSEA_"),
         starts_with("quanTIseq_"), starts_with("MCP_"),
         starts_with("EPIC_"), starts_with("xCell_"))
cat("Dimension of immune features for clustering:", dim(immune_features_for_clustering), "\n")#1013 125 

#Building a matrix for clustering analysis
clustering_matrix <- immune_features_for_clustering %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()
head(clustering_matrix[, 1:10])

# Calculate variance for each feature and select top 100 high-variance features

feature_variance <- apply(clustering_matrix, 2, var, na.rm = TRUE)
head(feature_variance)

high_var_features <- names(sort(feature_variance, decreasing = TRUE)[1:100])

clustering_features <- high_var_features

# count per method
cat("Feature breakdown:\n")
table(gsub("_.*", "", clustering_features)) 

## ---------------------------------------------------------------------------
## 4. Funtion for consensus clustering analysis 
## ---------------------------------------------------------------------------

run_consensus_clustering <- function(data, subtype_name, features, output_dir) {

  cat("\n=== Processing", subtype_name, "===\n")

  subtype_data <- data %>% filter(clinical_subtype == subtype_name)

  if (nrow(subtype_data) < 20) {
    cat("Warning: Too few samples (", nrow(subtype_data), ") for", subtype_name, "\n")
    return(NULL)
  }

  cat("Samples:", nrow(subtype_data), "\n")

  mat <- as.matrix(subtype_data[, features])
  rownames(mat) <- subtype_data$sample
  mat <- mat[, apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)]
  mat_scaled <- t(scale(t(mat)))

  results_dir <- file.path(output_dir, subtype_name)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  cc_results <- ConsensusClusterPlus(
    d         = mat_scaled,
    maxK      = 6,
    reps      = 1000,
    pItem     = 0.8,
    pFeature  = 0.8,
    clusterAlg = "hc",
    distance  = "euclidean",
    title     = results_dir,
    plot      = "png",
    verbose   = TRUE
  )

  return(list(
    consensus = cc_results,
    data      = subtype_data,
    matrix    = mat_scaled
  ))
}