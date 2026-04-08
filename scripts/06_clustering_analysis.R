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

# Merge immune features with clinical
data_main <- immune_features %>%
  left_join(clinical, by = c("sample" = "barcode")) %>%
  # Remove CIBERSORT QC columns (Used for validation later)
  select(-any_of(c("CIBERSORT_P-value",
                   "CIBERSORT_Correlation",
                   "CIBERSORT_RMSE")))

cat("Total samples loaded:", dim(data_main), "\n")
cat("PAM50 subtype distribution:\n")
print(table(data_main$PAM50_Subtype))

#===============================================================================
# 2. FEATURE SELECTION
# Use only RNA-seq validated deconvolution methods
# Exclude CIBERSORT from clustering (used later for validation)
# Select top 100 high-variance features to focus on informative signal
#===============================================================================

# Select features from validated methods only
immune_cols <- grep(
  "^(GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_)",
  colnames(data_main),
  value = TRUE
)

cat("\nTotal immune features available:", length(immune_cols), "\n")

# Build feature matrix: samples * features
immune_matrix <- as.matrix(data_main[, immune_cols])
rownames(immune_matrix) <- data_main$sample

# Remove zero-variance features (carry no information)
immune_matrix <- immune_matrix[,
  apply(immune_matrix, 2, function(x) var(x, na.rm = TRUE) > 0)
]

# Select top 100 most variable features across all samples

feature_variance    <- apply(immune_matrix, 2, var, na.rm = TRUE)
high_var_features   <- names(sort(feature_variance, decreasing = TRUE)[1:100])
clustering_features <- high_var_features

cat("Features selected (top variance):", length(clustering_features), "\n")
cat("Feature method breakdown:\n")
print(table(gsub("_.*", "", clustering_features)))

#===============================================================================
# 3. PCA REDUCTION FUNCTION
#===============================================================================

select_significant_pcs <- function(mat_scaled, subtype_name,
                                    var_threshold = 0.80) {

  # Run PCA on samples * features matrix
  pca_res <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  var_exp <- summary(pca_res)$importance[2, ]   # variance per PC
  cum_var <- summary(pca_res)$importance[3, ]   # cumulative variance

  # Select PCs that together explain var_threshold of total variance
  n_pcs <- which(cum_var >= var_threshold)[1]
  n_pcs <- max(n_pcs, 3)    # use at least 3 PCs
  n_pcs <- min(n_pcs, 20)   # cap at 20 PCs

  cat("  PCs selected:", n_pcs,
      "| Cumulative variance:", round(cum_var[n_pcs] * 100, 1), "%\n")

  # Scree plot — shows how much each PC contributes
  scree_df <- data.frame(
    PC      = seq_len(min(20, length(var_exp))),
    var_exp = var_exp[seq_len(min(20, length(var_exp)))] * 100,
    cum_var = cum_var[seq_len(min(20, length(var_exp)))] * 100
  )

  p_scree <- ggplot(scree_df, aes(x = PC)) +
    geom_bar(aes(y = var_exp), stat = "identity",
             fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = cum_var / 5), color = "red", linewidth = 1) +
    geom_point(aes(y = cum_var / 5), color = "red", size = 2) +
    geom_vline(xintercept = n_pcs, linetype = "dashed",
               color = "darkred", linewidth = 0.8) +
    scale_y_continuous(
      name     = "Variance Explained (%)",
      sec.axis = sec_axis(~ . * 5, name = "Cumulative Variance (%)")
    ) +
    scale_x_continuous(breaks = seq_len(20)) +
    labs(title    = paste(subtype_name, "— Scree Plot"),
         subtitle = paste(n_pcs, "PCs selected |",
                          round(cum_var[n_pcs] * 100, 1),
                          "% variance explained"),
         x = "Principal Component") +
    theme_bw()

  ggsave(
    paste0("results/figures/clustering/scree_", subtype_name, ".png"),
    p_scree, width = 8, height = 5, dpi = 300
  )

  return(list(
    pc_scores = pca_res$x[, seq_len(n_pcs), drop = FALSE],
    n_pcs     = n_pcs,
    var_exp   = var_exp,
    cum_var   = cum_var,
    pca_res   = pca_res
  ))
}

#===============================================================================
# 4. CONSENSUS CLUSTERING FUNCTION (PCA-based)
# ConsensusClusterPlus runs clustering 1000x on random subsamples
# This gives STABLE cluster assignments (not sensitive to random seed)
# pFeature=1.0 because we already reduced to PCs — use all of them
# d = t(pc_scores): function clusters COLUMNS, so columns must be samples
#===============================================================================

run_consensus_clustering_pca <- function(data, subtype_name,
                                          features, output_dir,
                                          var_threshold = 0.80) {

  cat("\n=== Processing", subtype_name, "===\n")

  subtype_data <- data %>% filter(PAM50_Subtype == subtype_name)

  if (nrow(subtype_data) < 20) {
    cat("Warning: Too few samples (", nrow(subtype_data), ") for",
        subtype_name, "\n")
    return(NULL)
  }

  cat("Samples:", nrow(subtype_data), "\n")

  # Build matrix: samples x features
  mat <- as.matrix(subtype_data[, features])
  rownames(mat) <- subtype_data$sample

  # Remove zero-variance features
  mat <- mat[, apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)]

  # Scale each feature across samples (z-score per feature)
  mat_scaled <- scale(mat)

  # Reduce dimensions with PCA
  cat("Running PCA...\n")
  pca_result <- select_significant_pcs(mat_scaled, subtype_name, var_threshold)
  pc_scores  <- pca_result$pc_scores   # samples x n_pcs

  cat("Clustering on", ncol(pc_scores), "PCs (reduced from",
      ncol(mat_scaled), "features)\n")

  results_dir <- file.path(output_dir, subtype_name)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # Run ConsensusClusterPlus
  # t(pc_scores) = PCs x samples (features x samples format required)
  cc_results <- ConsensusClusterPlus(
    d          = t(pc_scores),  # PCs x samples — clusters columns (samples)
    maxK       = 6,
    reps       = 1000,          # 1000 bootstrap iterations for stability
    pItem      = 0.8,           # 80% of samples per iteration
    pFeature   = 1.0,           # 100% of PCs (already reduced)
    clusterAlg = "hc",          # hierarchical clustering
    distance   = "euclidean",
    title      = results_dir,
    plot       = "png",
    verbose    = TRUE
  )

  # Verify clustering was done on SAMPLES not features
  n_consensus <- length(cc_results[[2]]$consensusClass)
  n_samples   <- nrow(pc_scores)
  cat("Verification — consensusClass length:", n_consensus,
      "| n_samples:", n_samples,
      "| Match:", n_consensus == n_samples, "\n")

  if (n_consensus != n_samples) {
    stop("Dimension mismatch in ", subtype_name,
         ". Check matrix orientation.")
  }

  return(list(
    consensus  = cc_results,
    data       = subtype_data,
    matrix     = mat_scaled,   # original scaled matrix for heatmaps
    pc_scores  = pc_scores,    # PC matrix for silhouette + UMAP
    pca_result = pca_result,
    n_pcs      = pca_result$n_pcs
  ))
}

#===============================================================================
# 5. OPTIMAL K DETERMINATION
# Why: We don't know how many immune subgroups exist within each subtype
#      Silhouette score measures how well-separated clusters are (-1 to +1)
#      Higher = better separation
#      We enforce minimum 10% in smallest cluster to avoid outlier detection
#===============================================================================

determine_optimal_k_pca <- function(cc_result, subtype_name,
                                     min_cluster_pct = 0.10) {

  cat("\n--- Optimal k for", subtype_name, "---\n")

  # Use PC scores for distance calculation
  # Why: Same space used for clustering → consistent distance metric
  mat       <- cc_result$pc_scores
  n_samples <- nrow(mat)
  dist_mat  <- dist(mat, method = "euclidean")

  cat("n_samples:", n_samples, "| n_PCs:", ncol(mat), "\n")

  k_range <- 2:6
  metrics  <- data.frame(
    k          = k_range,
    silhouette = NA_real_,
    sizes      = NA_character_,
    min_pct    = NA_real_,
    balanced   = NA
  )

  for (k in k_range) {

    clusters <- cc_result$consensus[[k]]$consensusClass

    # Guard checks — skip invalid results
    if (is.null(clusters) || any(is.na(clusters)))      { next }
    if (length(unique(clusters)) < 2)                   { next }
    if (length(clusters) != attr(dist_mat, "Size"))     { next }

    # Silhouette: how well does each sample fit its cluster?
    sil     <- silhouette(clusters, dist_mat)
    avg_sil <- mean(sil[, 3])

    # Cluster balance check
    tbl     <- table(clusters)
    sizes   <- paste(sort(as.numeric(tbl), decreasing = TRUE),
                     collapse = "/")
    min_pct <- min(as.numeric(tbl)) / n_samples * 100

    metrics$silhouette[metrics$k == k] <- avg_sil
    metrics$sizes[metrics$k == k]      <- sizes
    metrics$min_pct[metrics$k == k]    <- round(min_pct, 1)
    metrics$balanced[metrics$k == k]   <- min_pct >= (min_cluster_pct * 100)

    cat("  k =", k,
        "| Silhouette =", round(avg_sil, 3),
        "| Sizes:", sizes,
        "| Min%:", round(min_pct, 1), "%",
        "| Balanced:", min_pct >= (min_cluster_pct * 100), "\n")
  }

  # Select best k from BALANCED solutions only
  balanced_metrics <- metrics %>% filter(balanced == TRUE)

  if (nrow(balanced_metrics) == 0) {
    cat("  No balanced k found — data has no discrete immune subgroups\n")
    optimal_k    <- NA
    balance_flag <- "NO_SUBGROUPS"
  } else {
    optimal_k    <- balanced_metrics$k[which.max(balanced_metrics$silhouette)]
    balance_flag <- "OK"
  }

  cat("→ Optimal k:", optimal_k, "| Status:", balance_flag, "\n")

  # Plot silhouette by k
  p <- ggplot(metrics %>% filter(!is.na(silhouette)),
              aes(x = k, y = silhouette)) +
    geom_rect(
      data = metrics %>% filter(!is.na(balanced), !balanced),
      aes(xmin = k - 0.4, xmax = k + 0.4, ymin = -Inf, ymax = Inf),
      fill = "red", alpha = 0.1, inherit.aes = FALSE
    ) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(aes(color = balanced, shape = balanced), size = 4) +
    { if (!is.na(optimal_k))
        geom_vline(xintercept = optimal_k, color = "red",
                   linetype = "dashed", linewidth = 0.8) } +
    geom_text(aes(label = round(silhouette, 3)), vjust = -0.8, size = 3) +
    geom_text(aes(label = paste0(sizes, "\n(", min_pct, "%)")),
              vjust = 2.5, size = 2.5, color = "gray40") +
    scale_color_manual(
      values = c("TRUE" = "steelblue", "FALSE" = "red"),
      name   = "Balanced (>=10%)"
    ) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 17),
      name   = "Balanced (>=10%)"
    ) +
    scale_x_continuous(breaks = k_range) +
    labs(
      title    = paste(subtype_name, "— Optimal k Selection"),
      subtitle = paste("n =", n_samples, "| PCs:", ncol(mat),
                       "| Status:", balance_flag),
      x = "Number of Clusters (k)",
      y = "Mean Silhouette Width"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    paste0("results/figures/clustering/optimal_k_", subtype_name, ".png"),
    p, width = 8, height = 6, dpi = 300
  )

  return(list(
    metrics      = metrics,
    optimal_k    = optimal_k,
    balance_flag = balance_flag
  ))
}

#===============================================================================
# 6. RUN CLUSTERING FOR ALL SUBTYPES
#===============================================================================

subtypes <- c("Basal", "Her2", "LumA", "LumB", "Normal")

clustering_results_pca <- list()

for (subtype in subtypes) {
  clustering_results_pca[[subtype]] <- run_consensus_clustering_pca(
    data          = data_main,
    subtype_name  = subtype,
    features      = clustering_features,
    output_dir    = "results/consensus_clustering_pca",
    var_threshold = 0.80
  )
}

saveRDS(clustering_results_pca,
        "data/processed/clustering/consensus_clustering_pca.rds")

# Determine optimal k
optimal_k_pca <- list()

for (subtype in subtypes) {
  if (!is.null(clustering_results_pca[[subtype]])) {
    optimal_k_pca[[subtype]] <- determine_optimal_k_pca(
      clustering_results_pca[[subtype]],
      subtype,
      min_cluster_pct = 0.10
    )
  }
}

saveRDS(optimal_k_pca,
        "results/consensus_clustering_pca/optimal_k_pca.rds")



#===============================================================================
# 7. DECISION — DISCRETE CLUSTERS vs CONTINUOUS SCORES
#===============================================================================



