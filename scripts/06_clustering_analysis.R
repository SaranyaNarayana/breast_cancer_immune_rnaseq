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

table(data_full$PAM50_Subtype)


## -----------------------------
## 2. CIBERSORT  QC filtering
## -----------------------------
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

  data_full %>%
    select(`CIBERSORT_P-value`, `CIBERSORT_RMSE`,  `CIBERSORT_Correlation`) %>%
    summary()


# Apply QC filter
data_full_2 <- data_full %>%
  filter(`CIBERSORT_P-value` < 0.05,
         `CIBERSORT_Correlation` > 0,
         `CIBERSORT_RMSE` < quantile(`CIBERSORT_RMSE`, 0.95, na.rm = TRUE))
cat("After CIBERSORT QC filtering, retained", dim(data_full_2),"\n")

cat("PAM50 Subtype distribution before QC filtering:\n")
table(data_full$PAM50_Subtype)

cat("PAM50 Subtype distribution after QC filtering:\n")
table(data_full_2$PAM50_Subtype)

