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


# Remove zero-variance features
clustering_feature_cols <- clustering_feature_cols[
  apply(data_main[, clustering_feature_cols], 2,
        function(x) var(x, na.rm = TRUE) > 0)
]


log_message(paste("Clustering features:", length(clustering_feature_cols)))
log_message("Feature method breakdown:")
print(table(gsub("_.*", "", clustering_feature_cols)))

subtypes <- c("Basal", "Her2", "LumA", "LumB", "Normal")


## ---------------------------------------------------------------------------
## 5. Consensus Clustering within each PAM50 subtype
## ---------------------------------------------------------------------------


log_message("=== SECTION 2: CONSENSUS CLUSTERING ===")

# Global color palette — used across ALL figures for consistency
SUBTYPE_COLORS <- c(
  Basal  = "#E41A1C",
  Her2   = "#FF7F00",
  LumA   = "#4DAF4A",
  LumB   = "#377EB8",
  Normal = "#984EA3"
)

CLUSTER_COLORS <- c(
  IC1 = "#E41A1C",
  IC2 = "#377EB8",
  IC3 = "#4DAF4A",
  IC4 = "#FF7F00",
  IC5 = "#984EA3",
  IC6 = "#A65628"
)

GROUP_COLORS <- c("High" = "#E41A1C", "Low" = "#377EB8")


run_consensus_clustering <- function(data, subtype_name, features,
                                      output_dir, var_threshold = 0.80) {

  log_message(paste("\n--- Processing:", subtype_name, "---"))

  subtype_data <- data %>% filter(PAM50_Subtype == subtype_name)
  n_samples    <- nrow(subtype_data)
  log_message(paste("Samples:", n_samples))

  if (n_samples < 20) {
    log_message("WARNING: Too few samples (<20) — skipping")
    return(NULL)
  }

  # Build matrix: samples x features
  mat <- as.matrix(subtype_data[, features])
  rownames(mat) <- subtype_data$sample

  # Remove zero-variance features within this subtype
  mat <- mat[, apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)]
  log_message(paste("Features after variance filter:", ncol(mat)))

  # Z-score scale each feature across samples
  # WHY: Ensures no single feature dominates due to scale differences
  #      GSVA (-1 to 1) vs MCP (0 to 500) — without scaling MCP dominates
  mat_scaled <- scale(mat)

  # PCA reduction
  pca_res   <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  var_exp   <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cum_var   <- cumsum(var_exp)
  n_pcs     <- max(which(cum_var <= var_threshold), 3)
  n_pcs     <- min(n_pcs, 20)
  pc_scores <- pca_res$x[, seq_len(n_pcs), drop = FALSE]

  log_message(paste("PCs selected:", n_pcs,
                    "| Variance explained:", round(cum_var[n_pcs]*100, 1), "%"))

  # Save scree plot data
  scree_df <- data.frame(
    PC      = seq_len(min(20, length(var_exp))),
    var_exp = var_exp[seq_len(min(20, length(var_exp)))] * 100,
    cum_var = cum_var[seq_len(min(20, length(var_exp)))] * 100
  )

  p_scree <- ggplot(scree_df, aes(x = PC)) +
    geom_bar(aes(y = var_exp), stat = "identity",
             fill = SUBTYPE_COLORS[subtype_name], alpha = 0.7) +
    geom_line(aes(y = cum_var / 5), color = "black", linewidth = 1) +
    geom_point(aes(y = cum_var / 5), color = "black", size = 2) +
    geom_vline(xintercept = n_pcs, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_hline(yintercept = var_threshold * 100 / 5,
               linetype = "dotted", color = "gray50") +
    scale_y_continuous(
      name     = "Variance Explained per PC (%)",
      sec.axis = sec_axis(~ . * 5, name = "Cumulative Variance (%)")
    ) +
    scale_x_continuous(breaks = seq_len(20)) +
    labs(
      title    = paste(subtype_name, "— PCA Scree Plot"),
      subtitle = paste(n_pcs, "PCs selected |",
                       round(cum_var[n_pcs]*100, 1), "% total variance"),
      x = "Principal Component"
    ) +
    theme_bw()

  ggsave(
    paste0("results/figures/clustering/pca/scree_", subtype_name, ".png"),
    p_scree, width = 8, height = 5, dpi = 300
  )

  # Consensus clustering
  # t(pc_scores): CCP clusters COLUMNS → columns must be SAMPLES
  # pFeature = 1.0: use ALL PCs (already dimensionality-reduced)
  results_dir <- file.path(output_dir, subtype_name)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  cc_results <- ConsensusClusterPlus(
    d          = t(pc_scores),
    maxK       = 6,
    reps       = 1000,
    pItem      = 0.8,
    pFeature   = 1.0,
    clusterAlg = "hc",
    distance   = "euclidean",
    title      = results_dir,
    plot       = "png",
    verbose    = FALSE
  )

  # Verify orientation: consensusClass length must = n_samples
  n_consensus <- length(cc_results[[2]]$consensusClass)
  log_message(paste("Verification: consensusClass =", n_consensus,
                    "| n_samples =", n_samples,
                    "| Match =", n_consensus == n_samples))

  if (n_consensus != n_samples) {
    stop(paste("Dimension mismatch for", subtype_name,
               "— consensusClass:", n_consensus,
               "samples:", n_samples))
  }

  return(list(
    cc_results  = cc_results,
    data        = subtype_data,
    mat_scaled  = mat_scaled,
    pc_scores   = pc_scores,
    pca_res     = pca_res,
    var_exp     = var_exp,
    cum_var     = cum_var,
    n_pcs       = n_pcs,
    n_samples   = n_samples
  ))
}

clustering_results <- list()

for (subtype in subtypes) {
  clustering_results[[subtype]] <- run_consensus_clustering(
    data         = data_main,
    subtype_name = subtype,
    features     = clustering_feature_cols,
    output_dir   = "results/consensus_clustering",
    var_threshold = 0.80
  )
}

saveRDS(clustering_results,
        "data/processed/clustering/consensus_clustering_results.rds")
