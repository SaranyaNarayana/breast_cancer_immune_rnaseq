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
library(ConsensusClusterPlus)  # consensus clustering
library(cluster)               # silhouette scores
library(factoextra)            # PCA visualization
library(umap)                  # UMAP dimensionality reduction
library(dplyr)                 # data manipulation
library(ggplot2)               # plotting
library(pheatmap)              # heatmaps
library(RColorBrewer)          # color palettes
library(cowplot)               # combine plots
library(tidyr)                 # data reshaping
library(ggpubr)                # stat comparisons on plots
})


# Create output directories
dir.create("results/figures/clustering",        recursive = TRUE, showWarnings = FALSE)
dir.create("results/consensus_clustering_pca",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/clustering",                    recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/clustering",         recursive = TRUE, showWarnings = FALSE)

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
## 3. Define immune groups
## ---------------------------------------------------------------------------

#Remove ALL CIBERSORT columns from main data
data_main <- data_full %>%
  select(-starts_with("CIBERSORT_"))
cat("Dimension of data main after removing CIBERSORT columns:", dim(data_main), "\n")#1013 153


log_message("Defining immune component groups")
immune_component_groups <- list(

  # CD8 cytotoxic T cells + killing capacity
  # WHY: Core anti-tumor immune effectors
  #      Best predictor of immunotherapy response
  Cytotoxic = c(
    "GSVA_T_cells_CD8", "GSVA_Cytotoxicity",
    "ssGSEA_T_cells_CD8", "ssGSEA_Cytotoxicity",
    "quanTIseq_T.cells.CD8", "MCP_CD8 T cells", "EPIC_CD8_Tcells"
  ),

  # IFN signaling + NF-kB + TNF pathways
  # WHY: Broad inflammatory activation
  #      Includes your Hallmark MSigDB signatures from Script 05
  Inflammatory = c(
    "GSVA_HALLMARK_INFLAMMATORY_RESPONSE",
    "GSVA_HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "GSVA_HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "GSVA_HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "ssGSEA_HALLMARK_INFLAMMATORY_RESPONSE",
    "ssGSEA_HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "ssGSEA_HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "ssGSEA_HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  ),

  # Tregs + exhausted T cells + immunosuppressive pathways
  # WHY: Counter-regulatory mechanisms that shut down anti-tumor immunity
  #      High suppression = tumor is evading immune attack
  Suppression = c(
    "GSVA_T_reg", "GSVA_T_exhaustion", "GSVA_Immunosuppression",
    "ssGSEA_T_reg", "ssGSEA_T_exhaustion", "ssGSEA_Immunosuppression",
    "quanTIseq_Tregs"
  ),

  # M1 (classically activated, pro-inflammatory) macrophages
  # WHY: Anti-tumor macrophage phenotype
  #      Separate from M2 because they have OPPOSITE effects on survival
  Macrophage_M1 = c(
    "GSVA_Macrophages_M1", "ssGSEA_Macrophages_M1",
    "quanTIseq_Macrophages.M1"
  ),

  # M2 (alternatively activated, immunosuppressive) macrophages
  # WHY: Pro-tumor macrophage phenotype — promotes metastasis
  #      Kept separate from M1 — they are biologically opposite
  Macrophage_M2 = c(
    "GSVA_Macrophages_M2", "ssGSEA_Macrophages_M2",
    "quanTIseq_Macrophages.M2", "EPIC_Macrophages"
  ),

  # B cells and humoral immunity
  # WHY: B cells form tertiary lymphoid structures in some breast cancers
  #      Associated with good prognosis in TNBC
  B_cells = c(
    "GSVA_B_cells", "ssGSEA_B_cells", "quanTIseq_B.cells",
    "MCP_B lineage", "EPIC_Bcells"
  ),

  # Natural killer cells — innate cytotoxic lymphocytes
  # WHY: Innate immune surveillance — independent of antigen presentation
  NK_cells = c(
    "quanTIseq_NK.cells", "MCP_NK cells", "EPIC_NKcells"
  ),

  # Dendritic cells — antigen presenting cells
  # WHY: Required to prime T cells — gateway between innate and adaptive immunity
  Dendritic = c(
    "GSVA_Dendritic_cells", "ssGSEA_Dendritic_cells",
    "quanTIseq_Dendritic.cells", "MCP_Myeloid dendritic cells"
  ),

  # CAFs + endothelial — structural tumor microenvironment
  # WHY: Stromal cells physically exclude immune cells
  #      High stromal content = immune-excluded phenotype
  #      Independent from immune activation/suppression axis
  Stromal = c(
    "MCP_Fibroblasts", "MCP_Endothelial cells",
    "EPIC_CAFs", "EPIC_Endothelial"
  )
)

# Keep only features present in your data
immune_component_groups <- lapply(immune_component_groups, function(feats) {
  feats[feats %in% colnames(data_main)]
})
immune_component_groups <- immune_component_groups[
  sapply(immune_component_groups, length) > 0
]

log_message("Immune component groups:")
for (g in names(immune_component_groups)) {
  log_message(sprintf("  %-20s: %d features", g, length(immune_component_groups[[g]])))
}
