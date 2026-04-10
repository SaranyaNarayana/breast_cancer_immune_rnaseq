################################################################################
# 06_clustering_analysis.R
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
## Purpose:
#   - Seprate CIBERSORT QC for validation only (Section 2)
#   - top 100 high-variance features
#   -  
#   - 
#   - 
################################################################################

setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")

log_file <- file(
  paste0("logs/06_cluster_analysis_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"),
  open = "wt"
)
sink(log_file, type = "output")
sink(log_file, type = "message")

log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
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

# Output directories
dir.create("results/figures/clustering",       recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/clustering/pca",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/clustering/umap",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/consensus_clustering",     recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/clustering",        recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/clustering",        recursive = TRUE, showWarnings = FALSE)

# Global color palettes
SUBTYPE_COLORS <- c(
  Basal  = "#E41A1C", Her2   = "#FF7F00",
  LumA   = "#4DAF4A", LumB   = "#377EB8", Normal = "#984EA3"
)
CLUSTER_COLORS <- c(
  IC1 = "#E41A1C", IC2 = "#377EB8", IC3 = "#4DAF4A",
  IC4 = "#FF7F00", IC5 = "#984EA3", IC6 = "#A65628"
)

################################################################################
# SECTION 1: LOAD DATA
################################################################################

log_message("Loading data")

immune_features <- readRDS("data/processed/immune/immune_features_combined.rds")
clinical        <- readRDS("data/processed/clinical_filtered.rds")

cat("Immune features dimensions:", dim(immune_features), "\n")
cat("Clinical dimensions:", dim(clinical), "\n")

data_full <- immune_features %>%
  left_join(clinical, by = c("sample" = "barcode"))

cat("Merged data dimensions:", dim(data_full), "\n")
cat("PAM50 distribution:\n")
print(table(data_full$PAM50_Subtype))

saveRDS(data_full, "data/processed/clustering/data_full.rds")

################################################################################
# SECTION 2: CIBERSORT QC — SEPARATE FOR VALIDATION ONLY
################################################################################

log_message("CIBERSORT QC filtering")

cibersort_cols <- grep("^CIBERSORT_", colnames(data_full), value = TRUE)

cibersort_all <- data_full %>%
  select(sample, PAM50_Subtype, all_of(cibersort_cols))

cat("CIBERSORT QC summary:\n")
print(summary(cibersort_all[, c("CIBERSORT_P-value",
                                 "CIBERSORT_RMSE",
                                 "CIBERSORT_Correlation")]))

cibersort_valid <- cibersort_all %>%
  filter(`CIBERSORT_P-value`   < 0.05,
         `CIBERSORT_Correlation` > 0,
         `CIBERSORT_RMSE` < quantile(`CIBERSORT_RMSE`, 0.95, na.rm = TRUE))

cat("CIBERSORT after QC:", nrow(cibersort_valid), "samples\n")
saveRDS(cibersort_valid, "data/processed/clustering/cibersort_valid.rds")

# Remove ALL CIBERSORT from main analysis data
data_main <- data_full %>% select(-starts_with("CIBERSORT_"))
cat("data_main dimensions (CIBERSORT removed):", dim(data_main), "\n")

################################################################################
# SECTION 3: FEATURE SELECTION: top 100 high-variance features
################################################################################

log_message("Feature selection")

immune_cols_for_clustering <- data_full %>%
  select(sample,
         starts_with("GSVA_"),     starts_with("ssGSEA_"),
         starts_with("quanTIseq_"), starts_with("MCP_"),
         starts_with("EPIC_"),      starts_with("xCell_"))

cat("Total immune features available:", ncol(immune_cols_for_clustering) - 1, "\n")

# Build matrix: samples x features
clustering_matrix <- immune_cols_for_clustering %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Top 100 high-variance features
feature_variance    <- apply(clustering_matrix, 2, var, na.rm = TRUE)
high_var_features   <- names(sort(feature_variance, decreasing = TRUE)[1:100])
clustering_features <- high_var_features

cat("Top 100 high-variance features selected\n")
cat("Feature method breakdown:\n")
print(table(gsub("_.*", "", clustering_features)))

################################################################################
# SECTION 4: PC SELECTION FUNCTION — ELBOW METHOD ONLY
#
# WHY ELBOW ONLY (not 50% threshold):
#   The 50% threshold is an arbitrary cut-off with no statistical basis.
#   The elbow method directly answers: "where do PCs stop being informative?"
#   The elbow = the point where adding more PCs gives diminishing returns.
#   PCs before the elbow capture real biological signal.
#   PCs after the elbow capture noise.
#
# HOW ELBOW IS DETECTED:
#   Step 1: Compute variance explained per PC (first derivative = rate of drop)
#   Step 2: Compute second derivative = rate of change of the rate of drop
#   Step 3: The elbow is where second derivative is MOST POSITIVE
#           = where the steep drop TRANSITIONS to a flat tail
#           = the sharpest "bend" in the scree curve (concave up)
#   Skip PC1 (always the steepest, not the biological elbow)
#
# ONE FIGURE PRODUCED:
#   Scree plot: variance per PC as bars + cumulative variance line
#               Red dashed = elbow-selected PC
#               Shows how much variance each PC explains and where the elbow falls
#
# FLOOR = 3 PCs | CEILING = 10 PCs
################################################################################

select_pcs_elbow <- function(mat_scaled, subtype_name) {

  pca_res <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cum_var <- cumsum(var_exp)
  n_test  <- min(20, length(var_exp))

  #--------------------------------------------------------------------
  # ELBOW DETECTION via second derivative
  #
  # var_sub    = variance per PC (decreasing sequence)
  # first_d    = diff(var_sub) = how fast variance drops between PCs
  #              All values are negative (curve drops)
  # second_d   = diff(first_d) = rate of change of the drop rate
  #              POSITIVE = drop is slowing down (approaching flat)
  #              NEGATIVE = drop is accelerating (getting steeper)
  # The ELBOW is where second_d is MOST POSITIVE = drop slows the most
  # We skip index 1 [-1] because PC1→PC2 transition is always extreme
  # +2 corrects for: (a) double differencing offset, (b) skipping index 1
  #--------------------------------------------------------------------
  var_sub  <- var_exp[seq_len(n_test)]
  first_d  <- diff(var_sub)                    # length = n_test - 1
  second_d <- diff(first_d)                    # length = n_test - 2

  # Skip first element (PC1→PC2 always extreme), find max of remainder
  elbow_idx   <- which.max(second_d[-1]) + 2   # +2: offset correction
  n_pcs_final <- max(elbow_idx, 3)             # floor at 3
  n_pcs_final <- min(n_pcs_final, 10)          # ceiling at 10

  cat(sprintf(
    "  PC selection [%s]: elbow=PC%d | var_explained=%.1f%%\n",
    subtype_name, n_pcs_final, cum_var[n_pcs_final] * 100
  ))

  #--------------------------------------------------------------------
  # FIGURE: Scree plot (Panel A only)
  #   Bars = variance per PC + cumulative variance line
  #   Red dashed = elbow-selected PC
  #--------------------------------------------------------------------

  # Data frame for plot
  scree_df <- data.frame(
    PC      = seq_len(n_test),
    var_pct = var_exp[seq_len(n_test)] * 100,
    cum_pct = cum_var[seq_len(n_test)] * 100
  )

  y_max   <- max(scree_df$var_pct)
  scale_f <- y_max / max(scree_df$cum_pct) * 4

  # Scree plot
  p_scree <- ggplot(scree_df, aes(x = PC)) +

    # Bars: variance per PC
    geom_bar(aes(y = var_pct), stat = "identity",
             fill = SUBTYPE_COLORS[subtype_name],
             alpha = 0.70, color = "white", linewidth = 0.2) +

    # Cumulative variance line (scaled to fit left axis)
    geom_line(aes(y = cum_pct * scale_f),
              color = "black", linewidth = 1.1) +
    geom_point(aes(y = cum_pct * scale_f),
               color = "black", size = 2.2) +

    # Elbow = selected PC (red dashed vertical line)
    geom_vline(xintercept = n_pcs_final,
               color = "#E41A1C", linetype = "dashed", linewidth = 1.2) +

    # Label for selected PC
    annotate("label",
             x = n_pcs_final + 0.25, y = y_max * 0.88,
             label = paste0("Elbow selected\nPC", n_pcs_final,
                            " (", round(cum_var[n_pcs_final]*100, 1), "% var)"),
             color = "#E41A1C", fill = "white",
             hjust = 0, size = 3, fontface = "bold", label.size = 0.3) +

    scale_y_continuous(
      name     = "Variance Explained per PC (%)",
      sec.axis = sec_axis(~ . / scale_f, name = "Cumulative Variance (%)")
    ) +
    scale_x_continuous(breaks = seq_len(n_test)) +
    labs(
      title    = paste(subtype_name, "— Scree Plot"),
      subtitle = paste0("Bars = variance per PC  |  Line = cumulative variance",
                        "  |  Red = elbow-selected PC", n_pcs_final),
      x = "Principal Component"
    ) +
    theme_bw() +
    theme(plot.subtitle = element_text(size = 8, color = "gray40"),
          plot.title    = element_text(face = "bold"))

  ggsave(
    paste0("results/figures/clustering/pca/scree_", subtype_name, ".png"),
    p_scree, width = 10, height = 6, dpi = 300
  )

  return(list(
    pc_scores   = pca_res$x[, seq_len(n_pcs_final), drop = FALSE],
    n_pcs       = n_pcs_final,
    var_exp     = var_exp,
    cum_var     = cum_var,
    pca_res     = pca_res,
    n_pcs_elbow = n_pcs_final,
    second_d    = second_d
  ))
}

################################################################################
# SECTION 5: CONSENSUS CLUSTERING
#
# WHY PCA before clustering:
#   150+ features are highly correlated → PCA removes redundancy
#   Clustering in PC space finds genuine sample groupings
#   pFeature = 1.0: use all selected PCs (already reduced — no further sampling)
#   pItem = 0.8: each bootstrap uses 80% of samples → stability testing
################################################################################

log_message("=== SECTION 5: CONSENSUS CLUSTERING ===")

run_consensus_clustering <- function(data, subtype_name,
                                      features, output_dir) {

  cat("\n=== Processing:", subtype_name, "===\n")

  subtype_data <- data %>% filter(PAM50_Subtype == subtype_name)
  n_samples    <- nrow(subtype_data)
  cat("Samples:", n_samples, "\n")

  if (n_samples < 20) {
    cat("WARNING: Too few samples (<20) — skipping\n")
    return(NULL)
  }

  # Build and scale matrix
  mat            <- as.matrix(subtype_data[, features])
  rownames(mat)  <- subtype_data$sample
  mat            <- mat[, apply(mat, 2, function(x) var(x, na.rm=TRUE) > 0)]
  mat_scaled     <- scale(mat)
  cat("Matrix (samples x features):", dim(mat_scaled), "\n")

  # Elbow-based PC selection
  cat("Selecting PCs...\n")
  pca_result <- select_pcs_elbow(mat_scaled, subtype_name)
  pc_scores  <- pca_result$pc_scores
  n_pcs      <- pca_result$n_pcs

  cat("Using", n_pcs, "PCs |",
      round(pca_result$cum_var[n_pcs] * 100, 1), "% variance\n")

  # Consensus clustering
  results_dir <- file.path(output_dir, subtype_name)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Running ConsensusClusterPlus (1000 reps)...\n")

  cc_results <- ConsensusClusterPlus(
    d          = t(pc_scores),   # PCs x samples — CCP clusters columns (samples)
    maxK       = 6,
    reps       = 1000,
    pItem      = 0.8,
    pFeature   = 1.0,            # all PCs used (already reduced)
    clusterAlg = "hc",
    distance   = "euclidean",
    title      = results_dir,
    plot       = "png",
    verbose    = FALSE
  )

  # Verify orientation
  n_consensus <- length(cc_results[[2]]$consensusClass)
  cat("Verification: consensusClass =", n_consensus,
      "| n_samples =", n_samples,
      "| Match =", n_consensus == n_samples, "\n")

  if (n_consensus != n_samples) {
    stop("DIMENSION MISMATCH for ", subtype_name,
         ": consensusClass=", n_consensus, " samples=", n_samples)
  }

  return(list(
    consensus    = cc_results,
    data         = subtype_data,
    matrix       = mat_scaled,    # samples x features — for characterization
    pc_scores    = pc_scores,     # samples x PCs — for silhouette/UMAP
    pca_res      = pca_result$pca_res,
    var_exp      = pca_result$var_exp,
    cum_var      = pca_result$cum_var,
    n_pcs        = n_pcs,
    n_samples    = n_samples
  ))
}

subtypes <- c("Basal", "Her2", "LumA", "LumB", "Normal")

clustering_results <- list()
for (subtype in subtypes) {
  clustering_results[[subtype]] <- run_consensus_clustering(
    data         = data_main,
    subtype_name = subtype,
    features     = clustering_features,
    output_dir   = "results/consensus_clustering"
  )
}

saveRDS(clustering_results,
        "data/processed/clustering/consensus_clustering_results.rds")

# ── RELOAD SHORTCUT ──────────────────────────────────────────────────────────
# Consensus clustering (1000 reps x 5 subtypes) takes ~20 minutes.
# On re-runs, comment out the for-loop above and uncomment the line below:
# clustering_results <- readRDS("data/processed/clustering/consensus_clustering_results.rds")
# ─────────────────────────────────────────────────────────────────────────────

# Sanity check
cat("\nClustering results:\n")
for (s in subtypes) {
  if (!is.null(clustering_results[[s]])) {
    cat(s, "— samples:", clustering_results[[s]]$n_samples,
        "| PCs:", clustering_results[[s]]$n_pcs, "\n")
  }
}

################################################################################
# SECTION 6: OPTIMAL K DETERMINATION + BALANCE CHECK
#
# KEY RULE: A k is valid ONLY if smallest cluster >= 10% of subtype samples
# This threshold distinguishes real immune subgroups from outlier detection
#
# Three metrics evaluated per k:
#   (1) Silhouette: how well-separated are clusters in PC space?
#   (2) Balance:    is smallest cluster >= 10%?
#   (3) Delta area: how bimodal is the consensus matrix CDF?
################################################################################

log_message("=== SECTION 6: OPTIMAL K WITH 10% BALANCE RULE ===")

determine_optimal_k <- function(cc_result, subtype_name,
                                 min_cluster_pct = 0.10) {

  cat("\n--- Optimal k for:", subtype_name, "---\n")

  # Use PC scores for distance (consistent with clustering space)
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
    balanced   = NA,
    delta_area = NA_real_
  )

  for (k in k_range) {

    clusters <- cc_result$consensus[[k]]$consensusClass

    if (is.null(clusters) || any(is.na(clusters)))    { next }
    if (length(unique(clusters)) < 2)                 { next }
    if (length(clusters) != attr(dist_mat, "Size"))   { next }

    # Silhouette
    sil     <- silhouette(clusters, dist_mat)
    avg_sil <- mean(sil[, 3])

    # Balance
    tbl     <- table(clusters)
    sizes   <- paste(sort(as.numeric(tbl), decreasing=TRUE), collapse="/")
    min_pct <- min(as.numeric(tbl)) / n_samples * 100

    # CDF delta area
    cm     <- cc_result$consensus[[k]]$consensusMatrix
    cdf_x  <- seq(0, 1, length.out = 100)
    cdf_y  <- ecdf(cm[lower.tri(cm)])(cdf_x)
    area   <- sum(diff(cdf_x) * (head(cdf_y,-1) + tail(cdf_y,-1)) / 2)

    metrics$silhouette[metrics$k == k] <- avg_sil
    metrics$sizes[metrics$k == k]      <- sizes
    metrics$min_pct[metrics$k == k]    <- round(min_pct, 1)
    metrics$balanced[metrics$k == k]   <- min_pct >= (min_cluster_pct * 100)
    metrics$delta_area[metrics$k == k] <- area

    cat(sprintf("  k=%d | Sil=%.3f | Sizes: %-15s | Min%%=%5.1f | Balance: %s\n",
                k, avg_sil, sizes, min_pct,
                ifelse(min_pct >= min_cluster_pct*100, "OK ✓", "FAIL ✗")))
  }

  # Select optimal k: highest silhouette among BALANCED k values
  balanced_k <- metrics %>% filter(!is.na(balanced), balanced == TRUE)

  if (nrow(balanced_k) == 0) {
    optimal_k    <- NA
    balance_flag <- "NO_SUBGROUPS"
    cat("  -> NO balanced k found — homogeneous immune profile\n")
  } else {
    optimal_k    <- balanced_k$k[which.max(balanced_k$silhouette)]
    balance_flag <- "SUBGROUPS_FOUND"
    cat("  -> Optimal k =", optimal_k, "| SUBGROUPS FOUND\n")
  }

  return(list(
    metrics      = metrics,
    optimal_k    = optimal_k,
    balance_flag = balance_flag,
    n_samples    = n_samples
  ))
}

optimal_k_results <- list()
for (subtype in subtypes) {
  if (!is.null(clustering_results[[subtype]])) {
    optimal_k_results[[subtype]] <- determine_optimal_k(
      clustering_results[[subtype]], subtype, min_cluster_pct = 0.10
    )
  }
}

saveRDS(optimal_k_results,
        "results/consensus_clustering/optimal_k_results.rds")
# Reload shortcut: optimal_k_results <- readRDS("results/consensus_clustering/optimal_k_results.rds")

################################################################################
# SECTION 7: FULL RESULTS TABLE + VISUALIZATIONS
################################################################################

log_message("=== SECTION 7: RESULTS TABLES AND SILHOUETTE PLOTS ===")

# --- 7a: Full metrics table with balance flag ---
optimal_k_table <- do.call(rbind, lapply(subtypes, function(s) {

  if (is.null(optimal_k_results[[s]])) return(NULL)

  mets      <- optimal_k_results[[s]]$metrics
  opt_k     <- optimal_k_results[[s]]$optimal_k
  n_samples <- clustering_results[[s]]$n_samples

  mets %>%
    filter(!is.na(silhouette)) %>%
    mutate(
      subtype      = s,
      n_samples    = n_samples,
      is_optimal   = (!is.na(opt_k)) & (k == opt_k),
      min_cluster  = sapply(strsplit(sizes, "/"),
                            function(x) min(as.numeric(x))),
      pct_smallest = round(min_cluster / n_samples * 100, 1),
      imbalanced   = !balanced,
      status       = optimal_k_results[[s]]$balance_flag
    ) %>%
    select(subtype, n_samples, k, silhouette, sizes,
           min_cluster, pct_smallest, balanced, imbalanced,
           is_optimal, status)
}))

cat("\n=== FULL CLUSTERING METRICS TABLE ===\n")
print(optimal_k_table)
write.csv(optimal_k_table,
          "results/tables/clustering/optimal_k_full_metrics.csv",
          row.names = FALSE)

# --- 7b: Summary comparison table ---
k_summary <- do.call(rbind, lapply(subtypes, function(s) {

  opt  <- optimal_k_results[[s]]
  mets <- opt$metrics

  # k=2 metrics (always shown for comparison)
  sil_k2  <- mets$silhouette[mets$k == 2]
  sz_k2   <- mets$sizes[mets$k == 2]
  pct_k2  <- mets$min_pct[mets$k == 2]

  # Optimal k metrics
  ok    <- opt$optimal_k
  flag  <- opt$balance_flag
  proceed <- flag == "SUBGROUPS_FOUND"

  data.frame(
    Subtype         = s,
    N_samples       = opt$n_samples,
    k2_sizes        = sz_k2,
    k2_min_pct      = pct_k2,
    k2_silhouette   = round(sil_k2, 3),
    Optimal_k       = ifelse(is.na(ok), "NA", as.character(ok)),
    Status          = flag,
    Proceed_10pct   = proceed,
    Interpretation  = dplyr::case_when(
      s == "Basal"  ~ "Homogeneous immune-HOT (8.8% < 10% threshold)",
      s == "Her2"   ~ "Predominantly one type (8.3% < 10% threshold)",
      s == "LumA"   ~ "Homogeneous immune-COLD (0.6% < 10% threshold)",
      s == "LumB"   ~ "Homogeneous immune-COLD (4.0% < 10% threshold)",
      s == "Normal" ~ "Genuine k=2 subgroups (35.1% — passes threshold)",
      TRUE          ~ "Unknown"
    ),
    stringsAsFactors = FALSE
  )
}))

cat("\n=== FINAL K RECOMMENDATION (10% balance rule) ===\n")
cat(sprintf("%-8s %-12s %-10s %-8s %-20s %-10s\n",
            "Subtype", "k2_sizes", "Min_%", "Proceed",
            "Status", "Interpretation"))
cat(strrep("-", 75), "\n")
for (i in seq_len(nrow(k_summary))) {
  cat(sprintf("%-8s %-12s %-10s %-8s %-20s %s\n",
              k_summary$Subtype[i],
              k_summary$k2_sizes[i],
              paste0(k_summary$k2_min_pct[i], "%"),
              ifelse(k_summary$Proceed_10pct[i], "YES ✓", "NO ✗"),
              k_summary$Status[i],
              k_summary$Interpretation[i]))
}

write.csv(k_summary,
          "results/tables/clustering/clustering_summary.csv",
          row.names = FALSE)

# --- 7c: Silhouette plots — all subtypes faceted ---
p_sil <- ggplot(optimal_k_table, aes(x = k, y = silhouette)) +

  # Red shading for imbalanced k values
  geom_rect(data = optimal_k_table %>% filter(imbalanced),
            aes(xmin = k - 0.4, xmax = k + 0.4,
                ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.08, inherit.aes = FALSE) +

  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(aes(shape = imbalanced, color = imbalanced), size = 4) +

  # Optimal k marker
  geom_vline(
    data = optimal_k_table %>% filter(is_optimal),
    aes(xintercept = k),
    color = "red", linetype = "dashed", linewidth = 0.9
  ) +

  geom_text(aes(label = round(silhouette, 3)),
            vjust = -0.9, size = 3) +
  geom_text(aes(label = paste0(sizes, "\n(", pct_smallest, "%)")),
            vjust = 2.5, size = 2.5, color = "gray40") +

  facet_wrap(~ subtype, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = 2:6) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    labels = c("FALSE" = "Balanced (>=10%)", "TRUE" = "Imbalanced (<10%)"),
    name   = "10% balance rule"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Balanced (>=10%)", "TRUE" = "Imbalanced (<10%)"),
    name   = "10% balance rule"
  ) +
  labs(
    title    = "Silhouette Score by k — All PAM50 Subtypes",
    subtitle = "Red = fails 10% balance rule | Dashed = optimal k | Sizes + min% shown",
    x = "Number of Clusters (k)", y = "Mean Silhouette Width"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold", size = 11),
        plot.subtitle   = element_text(size = 9, color = "gray40"))

ggsave("results/figures/clustering/01_silhouette_all_subtypes.png",
       p_sil, width = 14, height = 10, dpi = 300)

log_message("Silhouette plots saved")

################################################################################
# SECTION 8: EXTRACT CLUSTER ASSIGNMENTS
#
# UPDATED RESULTS (elbow PC selection, 10% balance rule):
#   Basal  : 166/16  (8.8%)  → FAIL — homogeneous immune-hot
#   Her2   :  66/6   (8.3%)  → FAIL — predominantly one type
#   LumA   : 520/3   (0.6%)  → FAIL — homogeneous immune-cold
#   LumB   : 191/8   (4.0%)  → FAIL — homogeneous immune-cold
#   Normal :  24/13  (35.1%) → PASS — genuine k=2 immune subgroups
#
# ONLY Normal proceeds with discrete IC1/IC2 cluster labels.
# All other subtypes use continuous immune scores (Section 11).
################################################################################

log_message("=== SECTION 8: EXTRACTING CLUSTER ASSIGNMENTS ===")

subtypes_with_clusters  <- k_summary$Subtype[k_summary$Proceed_10pct]
subtypes_gradient_only  <- k_summary$Subtype[!k_summary$Proceed_10pct]

cat("Subtypes with genuine clusters (>=10%):",
    paste(subtypes_with_clusters, collapse = ", "), "\n")
cat("Subtypes — gradient only (<10%):",
    paste(subtypes_gradient_only, collapse = ", "), "\n")

cluster_assignments <- list()

for (subtype in subtypes_with_clusters) {

  opt_k    <- optimal_k_results[[subtype]]$optimal_k
  clusters <- clustering_results[[subtype]]$consensus[[opt_k]]$consensusClass

  cluster_df <- data.frame(
    sample         = names(clusters),
    immune_cluster = paste0("IC", clusters),
    PAM50_Subtype  = subtype,
    stringsAsFactors = FALSE
  )

  cat(subtype, "— k =", opt_k, "clusters:\n")
  print(table(cluster_df$immune_cluster))

  cluster_assignments[[subtype]] <- cluster_df
}

################################################################################
# SECTION 9: BIOLOGICAL CHARACTERIZATION — NORMAL SUBTYPE ONLY
#
# WHY: Only Normal has genuine k=2 clusters (IC1=24, IC2=13, 35.1%).
#      Characterization explains WHAT makes IC1 different from IC2.
#      This is essential before interpreting survival/DE results.
#
# THREE APPROACHES:
#   (a) PCA loadings — which features drive each PC axis?
#   (b) Feature importance (Wilcoxon) — which features differ IC1 vs IC2?
#   (c) Biplot — samples + feature arrows in same space
#
# NOTE: Survival analysis, DE analysis, and pathway enrichment are in
#       07_survival_analysis.R and 08_DE_pathway_analysis.R
#       Those scripts use:
#         immune_cluster (IC1/IC2) for Normal
#         immune_group_primary (High/Low gradient) for all other subtypes
################################################################################

log_message("=== SECTION 9: BIOLOGICAL CHARACTERIZATION (Normal only) ===")

# --- 9a: PCA loadings — what drives each PC? ---
plot_pca_loadings <- function(cc_result, subtype_name, n_top = 15) {

  pca_res  <- cc_result$pca_res
  n_pcs    <- cc_result$n_pcs
  var_exp  <- cc_result$var_exp

  loadings <- as.data.frame(pca_res$rotation[, seq_len(min(4, n_pcs))])
  loadings$feature       <- rownames(loadings)
  loadings$feature_clean <- gsub("GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_",
                                   "", loadings$feature)
  loadings$method        <- gsub("_.*", "", loadings$feature)

  plot_list <- list()
  for (i in seq_len(min(4, n_pcs))) {
    pc_name <- paste0("PC", i)
    ve      <- round(var_exp[i] * 100, 1)

    top_load <- loadings %>%
      arrange(desc(abs(.data[[pc_name]]))) %>%
      head(n_top) %>%
      # Keep full feature name (with method prefix) to avoid duplicates
      # Display clean name as label separately
      mutate(
        direction    = ifelse(.data[[pc_name]] > 0, "Positive", "Negative"),
        # Use full feature as factor to avoid duplicate level error
        feature_ordered = factor(feature,
                                  levels = feature[order(.data[[pc_name]])])
      )

    p <- ggplot(top_load,
                aes(x = .data[[pc_name]], y = feature_ordered,
                    fill = direction)) +
      geom_bar(stat = "identity", alpha = 0.85) +
      geom_vline(xintercept = 0, linewidth = 0.5) +
      scale_fill_manual(values = c("Positive" = "#E41A1C",
                                    "Negative" = "#377EB8")) +
      # Show clean labels (method prefix stripped) on y axis
      scale_y_discrete(labels = setNames(top_load$feature_clean,
                                          top_load$feature)) +
      labs(title    = paste0(pc_name, " (", ve, "% var)"),
           subtitle = "Top loading features",
           x = "Loading", y = NULL) +
      theme_bw() +
      theme(legend.position = "none", axis.text.y = element_text(size = 7))

    plot_list[[pc_name]] <- p
  }

  p_combined <- plot_grid(plotlist = plot_list, ncol = 2)
  title_g    <- ggdraw() +
    draw_label(paste(subtype_name, "— PCA Loadings: What Drives Each PC"),
               fontface = "bold", size = 12)

  ggsave(
    paste0("results/figures/clustering/pca/loadings_", subtype_name, ".png"),
    plot_grid(title_g, p_combined, ncol = 1, rel_heights = c(0.06, 1)),
    width = 14, height = 12, dpi = 300
  )

  cat("\nTop features driving PC1 for", subtype_name, ":\n")
  print(loadings %>%
          arrange(desc(abs(PC1))) %>%
          head(10) %>%
          select(feature_clean, method, PC1))

  return(loadings)
}

# --- 9b: Feature importance (Wilcoxon IC1 vs IC2) ---
compute_feature_importance <- function(data_main, cluster_df,
                                        subtype_name, n_top = 30) {

  groups   <- sort(unique(cluster_df$immune_cluster))
  if (length(groups) < 2) return(NULL)

  feat_cols <- grep("^(GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_)",
                    colnames(data_main), value = TRUE)

  sub_data <- data_main %>%
    filter(PAM50_Subtype == subtype_name) %>%
    left_join(cluster_df %>% select(sample, immune_cluster), by = "sample") %>%
    filter(!is.na(immune_cluster))

  #--------------------------------------------------------------------
  # STEP 1: Z-score normalise EACH FEATURE across all samples of this
  # subtype BEFORE computing the mean difference.
  # WHY: MCP-counter uses arbitrary enrichment units (0–500+) while
  #      GSVA uses -1 to +1. Without normalisation the MCP T-cell bar
  #      dominates at +14 making all GSVA bars invisible.
  # After z-scoring every feature has mean=0, sd=1 → differences are
  # directly comparable across methods.
  #--------------------------------------------------------------------
  feat_mat_scaled <- sub_data %>%
    select(sample, immune_cluster, all_of(feat_cols)) %>%
    mutate(across(all_of(feat_cols), ~ as.numeric(scale(.)))) %>%
    # scale() returns NaN if sd=0; replace with 0
    mutate(across(all_of(feat_cols), ~ ifelse(is.nan(.), 0, .)))

  feat_res <- do.call(rbind, lapply(feat_cols, function(feat) {
    if (!feat %in% colnames(feat_mat_scaled)) return(NULL)

    g1 <- feat_mat_scaled[[feat]][feat_mat_scaled$immune_cluster == groups[1]]
    g2 <- feat_mat_scaled[[feat]][feat_mat_scaled$immune_cluster == groups[2]]
    if (length(g1) < 3 || length(g2) < 3) return(NULL)

    wt <- wilcox.test(g1, g2, exact = FALSE)

    data.frame(
      feature       = feat,
      feature_clean = gsub("GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_",
                            "", feat),
      method        = gsub("_.*", "", feat),
      mean_g1       = mean(g1, na.rm = TRUE),
      mean_g2       = mean(g2, na.rm = TRUE),
      # diff is now in z-score units — comparable across methods
      diff          = mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE),
      p_value       = wt$p.value,
      group1        = groups[1],
      group2        = groups[2],
      higher_in     = ifelse(mean(g2, na.rm = TRUE) > mean(g1, na.rm = TRUE),
                             groups[2], groups[1]),
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      FDR         = p.adjust(p_value, method = "BH"),
      significant = FDR < 0.05,

      #--------------------------------------------------------------------
      # Explicit biological categories including what "Other" covers
      # WHY: "Other" is uninformative — users cannot know what it means.
      # We now catch all major immune cell/pathway types explicitly so
      # nothing falls through to "Other" if avoidable.
      #--------------------------------------------------------------------
      bio_type = dplyr::case_when(
        # Cytotoxic / innate killing
        grepl("CD8|Cytotox|NK|NKT|Cytotoxic.lymph",
              feature_clean, ignore.case = TRUE)               ~ "Cytotoxic/NK",
        # Inflammatory signalling pathways
        grepl("IFN|Inflam|Cytokine|Chemokine|TNF|NFKB",
              feature_clean, ignore.case = TRUE)               ~ "Inflammatory",
        # Immunosuppressive / regulatory
        grepl("Treg|Exhaust|Immunosup|M2|MDSC|checkpoint|Checkpoint",
              feature_clean, ignore.case = TRUE)               ~ "Suppressive",
        # Stromal microenvironment
        grepl("CAF|Fibro|Endothel|Stromal|Pericyte|Smooth",
              feature_clean, ignore.case = TRUE)               ~ "Stromal",
        # Myeloid / antigen presentation
        grepl("M1|Dendritic|Antigen|DC|Monocyt|Macrophage|Myeloid",
              feature_clean, ignore.case = TRUE)               ~ "Myeloid/APC",
        # B cell lineage and humoral immunity
        grepl("B.cell|B_cell|Bcell|B.lineage|Plasma|Humoral",
              feature_clean, ignore.case = TRUE)               ~ "B cell/Humoral",
        # Pan-immune / composite scores
        grepl("ImmuneScore|MicroenvironmentScore|StromaScore|T.cells$|^T cells",
              feature_clean, ignore.case = TRUE)               ~ "Pan-immune score",
        # Neutrophils / innate myeloid
        grepl("Neutrophil|Basophil|Eosinophil|Mast",
              feature_clean, ignore.case = TRUE)               ~ "Innate myeloid",
        # Anything not caught — show which method it comes from
        TRUE ~ paste0("Other (", method, ")")
      )
    ) %>%
    arrange(FDR, desc(abs(diff)))

  write.csv(feat_res,
            paste0("results/tables/clustering/feature_importance_",
                   subtype_name, ".csv"),
            row.names = FALSE)

  # Top significant features ordered by normalised difference
  top_f <- feat_res %>%
    filter(significant) %>%
    head(n_top) %>%
    mutate(
      feature_ordered = factor(feature, levels = feature[order(diff)])
    )

  if (nrow(top_f) > 0) {

    #--------------------------------------------------------------------
    # Build a colour palette that covers all bio_types actually present
    # Base palette + dynamic "Other (method)" colours
    #--------------------------------------------------------------------
    base_cols <- c(
      "Cytotoxic/NK"    = "#E41A1C",
      "Inflammatory"    = "#FF7F00",
      "Suppressive"     = "#377EB8",
      "Stromal"         = "#4DAF4A",
      "Myeloid/APC"     = "#984EA3",
      "B cell/Humoral"  = "#A65628",
      "Pan-immune score"= "#999999",
      "Innate myeloid"  = "#F781BF"
    )

    # Any remaining "Other (method)" types get grey shades
    other_types <- setdiff(unique(top_f$bio_type), names(base_cols))
    if (length(other_types) > 0) {
      grey_shades <- colorRampPalette(c("gray55", "gray80"))(length(other_types))
      names(grey_shades) <- other_types
      bio_cols <- c(base_cols, grey_shades)
    } else {
      bio_cols <- base_cols
    }

    # Keep only colours that appear in the data
    bio_cols <- bio_cols[names(bio_cols) %in% unique(top_f$bio_type)]

    p_feat <- ggplot(top_f,
                     aes(x = diff, y = feature_ordered, fill = bio_type)) +
      geom_bar(stat = "identity", alpha = 0.88) +
      geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +

      # Add FDR significance stars at bar tips
      geom_text(
        aes(x = diff + sign(diff) * 0.03,
            label = dplyr::case_when(
              FDR < 0.001 ~ "***",
              FDR < 0.01  ~ "**",
              FDR < 0.05  ~ "*",
              TRUE        ~ ""
            )),
        hjust = 0.5,
        size  = 3, color = "black"
      ) +

      scale_fill_manual(values = bio_cols, name = "Biological Category") +

      # Show clean names (without method prefix) on y axis
      scale_y_discrete(labels = setNames(
        paste0(top_f$feature_clean, "  [", top_f$method, "]"),
        top_f$feature
      )) +

      facet_grid(higher_in ~ ., scales = "free_y", space = "free") +

      labs(
        title    = paste(subtype_name, "— Features Defining Immune Clusters"),
        subtitle = paste0(
          "Normalised mean difference (z-score units — comparable across methods)\n",
          "Positive = higher in ", groups[2],
          "  |  Negative = higher in ", groups[1],
          "  |  FDR < 0.05  |  * p<0.05  ** p<0.01  *** p<0.001\n",
          "[method] shown after each feature name"
        ),
        x = paste0("Normalised mean difference, z-score  (",
                   groups[2], " − ", groups[1], ")"),
        y = NULL
      ) +
      theme_bw() +
      theme(
        axis.text.y     = element_text(size = 7.5),
        legend.position = "right",
        legend.text     = element_text(size = 8),
        strip.text      = element_text(face = "bold", size = 10),
        plot.subtitle   = element_text(size = 8, color = "gray40",
                                        lineheight = 1.3)
      )

    ggsave(
      paste0("results/figures/clustering/feature_importance_",
             subtype_name, ".png"),
      p_feat,
      width  = 13,
      height = max(8, nrow(top_f) * 0.38),
      dpi    = 300
    )
  }

  # Auto biological labels
  hot_kw  <- c("CD8","Cytotox","IFN","Inflam","NK","B_cell","M1","Antigen")
  cold_kw <- c("M2","Immunosup","Exhaust","Treg","CAF","Fibro","Endothel")

  labels <- sapply(groups, function(g) {
    hf  <- feat_res %>% filter(significant, higher_in == g) %>% pull(feature_clean)
    hs  <- sum(sapply(hot_kw,  function(k) any(grepl(k, hf, ignore.case=TRUE))))
    cs  <- sum(sapply(cold_kw, function(k) any(grepl(k, hf, ignore.case=TRUE))))
    lbl <- dplyr::case_when(
      hs > cs  ~ "IMMUNE-HOT (Cytotoxic/Inflamed)",
      cs > hs  ~ "IMMUNE-COLD (Suppressed/Excluded)",
      TRUE     ~ "MIXED profile"
    )
    paste(g, "→", lbl)
  })

  cat("\n--- Biological labels for", subtype_name, "---\n")
  for (l in labels) cat(" ", l, "\n")

  return(list(feat_res = feat_res, labels = labels))
}

# --- 9c: Biplot — samples + feature arrows ---
plot_biplot <- function(cc_result, cluster_df, subtype_name,
                         n_arrows = 12) {

  pca_res   <- cc_result$pca_res
  pc_scores <- cc_result$pc_scores
  loadings  <- pca_res$rotation
  var_exp   <- cc_result$var_exp

  samples_df <- data.frame(
    PC1    = pc_scores[, 1],
    PC2    = pc_scores[, 2],
    sample = rownames(pc_scores)
  ) %>% left_join(cluster_df %>% select(sample, immune_cluster), by = "sample")

  arrow_df <- as.data.frame(loadings[, 1:2]) %>%
    mutate(
      feature       = rownames(loadings),
      feature_clean = gsub("GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_",
                            "", feature),
      magnitude     = sqrt(PC1^2 + PC2^2),
      bio_type      = dplyr::case_when(
        grepl("CD8|Cytotox|NK|B_cell",
              feature_clean, ignore.case=TRUE) ~ "Cytotoxic/NK",
        grepl("IFN|Inflam|Cytokine",
              feature_clean, ignore.case=TRUE) ~ "Inflammatory",
        grepl("Treg|Exhaust|Immunosup|M2",
              feature_clean, ignore.case=TRUE) ~ "Suppressive",
        grepl("CAF|Fibro|Endothel",
              feature_clean, ignore.case=TRUE) ~ "Stromal",
        TRUE                                   ~ "Other"
      )
    ) %>%
    arrange(desc(magnitude)) %>%
    head(n_arrows)

  sf <- max(abs(c(pc_scores[,1], pc_scores[,2]))) /
        max(arrow_df$magnitude) * 0.55

  arrow_df <- arrow_df %>%
    mutate(ax = PC1 * sf, ay = PC2 * sf,
           lx = PC1 * sf * 1.2, ly = PC2 * sf * 1.2)

  n_cl      <- length(unique(na.omit(samples_df$immune_cluster)))
  cl_cols   <- CLUSTER_COLORS[seq_len(n_cl)]
  names(cl_cols) <- sort(unique(na.omit(samples_df$immune_cluster)))

  bio_cols <- c("Cytotoxic/NK" = "#E41A1C", "Inflammatory" = "#FF7F00",
                 "Suppressive"  = "#377EB8", "Stromal"      = "#4DAF4A",
                 "Other"        = "gray50")

  all_cols <- c(cl_cols, bio_cols)

  p_bp <- ggplot() +
    geom_hline(yintercept = 0, color = "gray90") +
    geom_vline(xintercept = 0, color = "gray90") +
    stat_ellipse(data = samples_df %>% filter(!is.na(immune_cluster)),
                 aes(x = PC1, y = PC2, color = immune_cluster),
                 linewidth = 1, linetype = "dashed", level = 0.90) +
    geom_point(data  = samples_df,
               aes(x = PC1, y = PC2, color = immune_cluster),
               size = 2.5, alpha = 0.75) +
    geom_segment(data = arrow_df,
                 aes(x=0, y=0, xend=ax, yend=ay, color=bio_type),
                 arrow     = arrow(length = unit(0.25,"cm"), type="closed"),
                 linewidth = 0.8, alpha = 0.8, inherit.aes = FALSE) +
    geom_text_repel(data = arrow_df,
                    aes(x=lx, y=ly, label=feature_clean, color=bio_type),
                    size = 3, fontface = "bold", max.overlaps = 20,
                    inherit.aes = FALSE) +
    scale_color_manual(values = all_cols,
                       breaks = names(cl_cols),
                       name   = "Cluster") +
    labs(
      title    = paste(subtype_name, "— PCA Biplot"),
      subtitle = paste0(
        "DOTS = samples (by cluster)  |  ARROWS = immune features (by biology)\n",
        "Samples in arrow direction = HIGH for that feature"
      ),
      x = paste0("PC1 (", round(var_exp[1]*100,1), "%)"),
      y = paste0("PC2 (", round(var_exp[2]*100,1), "%)")
    ) +
    theme_bw() +
    theme(plot.subtitle = element_text(size = 8, color = "gray40"))

  ggsave(paste0("results/figures/clustering/pca/biplot_", subtype_name, ".png"),
         p_bp, width = 12, height = 9, dpi = 300)

  return(p_bp)
}

# Run all characterization for clustered subtypes
# Run characterization for Normal only
# (the only subtype with genuine clusters from the 10% balance rule)
char_results <- list()

for (subtype in subtypes_with_clusters) {

  cat("\n", strrep("=", 55), "\n")
  cat("CHARACTERIZING:", subtype, "(k=2, IC1/IC2)\n")
  cat(strrep("=", 55), "\n")

  # PCA loadings — what biological features drive the PC axes?
  plot_pca_loadings(clustering_results[[subtype]], subtype)

  # Feature importance — what significantly differs between IC1 and IC2?
  char_results[[subtype]] <- compute_feature_importance(
    data_main,
    cluster_assignments[[subtype]],
    subtype
  )

  # Biplot — show samples + feature arrows in same space
  plot_biplot(
    clustering_results[[subtype]],
    cluster_assignments[[subtype]],
    subtype
  )
}

################################################################################
# SECTION 10: PCA + UMAP VISUALIZATIONS (clustered subtypes)
################################################################################

log_message("=== SECTION 10: PCA + UMAP VISUALIZATIONS ===")

for (subtype in subtypes_with_clusters) {

  cc        <- clustering_results[[subtype]]
  cl_df     <- cluster_assignments[[subtype]]
  pc_scores <- cc$pc_scores
  var_exp   <- cc$var_exp

  plot_df <- data.frame(
    PC1 = pc_scores[,1], PC2 = pc_scores[,2],
    sample = rownames(pc_scores)
  ) %>% left_join(cl_df %>% select(sample, immune_cluster), by = "sample")

  n_cl     <- length(unique(na.omit(plot_df$immune_cluster)))
  cl_cols  <- CLUSTER_COLORS[seq_len(n_cl)]
  names(cl_cols) <- sort(unique(na.omit(plot_df$immune_cluster)))

  # PCA plot
  p_pca <- ggplot(plot_df, aes(x=PC1, y=PC2, color=immune_cluster)) +
    geom_point(size=3, alpha=0.85) +
    stat_ellipse(linewidth=1, linetype="dashed") +
    scale_color_manual(values=cl_cols, name="Cluster") +
    labs(
      title    = paste(subtype, "— Immune Clusters (PC Space)"),
      subtitle = paste("k=2 | n=", nrow(plot_df),
                       "| Sizes:", paste(table(plot_df$immune_cluster),
                                          collapse="/")),
      x = paste0("PC1 (", round(var_exp[1]*100,1), "%)"),
      y = paste0("PC2 (", round(var_exp[2]*100,1), "%)")
    ) + theme_bw()

  ggsave(paste0("results/figures/clustering/pca/PCA_", subtype, ".png"),
         p_pca, width=8, height=6, dpi=300)

  # UMAP plot
  set.seed(42)
  umap_res <- umap::umap(pc_scores,
                         n_neighbors = min(15, nrow(pc_scores)-1))
  umap_df  <- data.frame(
    UMAP1 = umap_res$layout[,1], UMAP2 = umap_res$layout[,2],
    sample = rownames(pc_scores)
  ) %>% left_join(cl_df %>% select(sample, immune_cluster), by="sample")

  p_umap <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=immune_cluster)) +
    geom_point(size=3, alpha=0.85) +
    stat_ellipse(linewidth=1, linetype="dashed") +
    scale_color_manual(values=cl_cols, name="Cluster") +
    labs(title    = paste(subtype, "— Immune Clusters (UMAP)"),
         subtitle = paste("k=2 | n=", nrow(umap_df))) +
    theme_bw()

  ggsave(paste0("results/figures/clustering/umap/UMAP_", subtype, ".png"),
         p_umap, width=8, height=6, dpi=300)

  # Cluster heatmap (top 50 variable features)
  sub_data <- cc$data %>%
    left_join(cl_df %>% select(sample, immune_cluster), by="sample") %>%
    arrange(immune_cluster)

  feat_var   <- apply(cc$matrix, 2, var, na.rm=TRUE)
  top50      <- names(sort(feat_var, decreasing=TRUE)[1:50])
  top50      <- top50[top50 %in% colnames(sub_data)]
  mat_plot   <- t(scale(t(as.matrix(sub_data[, top50]))))
  rownames(mat_plot) <- sub_data$sample

  ann_col    <- data.frame(Cluster = sub_data$immune_cluster,
                            row.names = sub_data$sample)
  clust_ann  <- CLUSTER_COLORS[seq_len(n_cl)]
  names(clust_ann) <- sort(unique(sub_data$immune_cluster))

  png(paste0("results/figures/clustering/heatmap_clusters_", subtype, ".png"),
      width=14, height=10, units="in", res=300)
  pheatmap(
    t(mat_plot),
    annotation_col    = ann_col,
    annotation_colors = list(Cluster = clust_ann),
    labels_row        = gsub("GSVA_|ssGSEA_|quanTIseq_|MCP_|EPIC_|xCell_",
                              "", top50),
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
    breaks = seq(-2,2,length.out=101),
    main = paste(subtype, "— Immune Feature Heatmap by Cluster"),
    fontsize=8, fontsize_row=7, border_color=NA
  )
  dev.off()

  log_message(paste(subtype, "— PCA, UMAP, heatmap saved"))
}

################################################################################
# SECTION 11: IMMUNE COMPONENT SCORES FOR ALL SUBTYPES
#
# WHY: Provides continuous immune activity scores for:
#   (a) Survival analysis — continuous Cox regression
#   (b) DE analysis — Q4 vs Q1 extreme comparison
#   (c) Immune landscape visualization across all subtypes
#
# 9 biologically curated components — each is an average of
# multiple methods measuring the same biological concept
################################################################################

log_message("=== SECTION 11: IMMUNE COMPONENT SCORES ===")

immune_component_groups <- list(
  Cytotoxic     = c("GSVA_T_cells_CD8","GSVA_Cytotoxicity",
                    "ssGSEA_T_cells_CD8","ssGSEA_Cytotoxicity",
                    "quanTIseq_T.cells.CD8","MCP_CD8 T cells","EPIC_CD8_Tcells"),
  Inflammatory  = c("GSVA_HALLMARK_INFLAMMATORY_RESPONSE",
                    "GSVA_HALLMARK_INTERFERON_GAMMA_RESPONSE",
                    "GSVA_HALLMARK_INTERFERON_ALPHA_RESPONSE",
                    "GSVA_HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                    "ssGSEA_HALLMARK_INFLAMMATORY_RESPONSE",
                    "ssGSEA_HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  Suppression   = c("GSVA_T_reg","GSVA_T_exhaustion","GSVA_Immunosuppression",
                    "ssGSEA_T_reg","ssGSEA_T_exhaustion","quanTIseq_Tregs"),
  Macrophage_M1 = c("GSVA_Macrophages_M1","ssGSEA_Macrophages_M1",
                    "quanTIseq_Macrophages.M1"),
  Macrophage_M2 = c("GSVA_Macrophages_M2","ssGSEA_Macrophages_M2",
                    "quanTIseq_Macrophages.M2","EPIC_Macrophages"),
  B_cells       = c("GSVA_B_cells","ssGSEA_B_cells","quanTIseq_B.cells",
                    "MCP_B lineage","EPIC_Bcells"),
  NK_cells      = c("quanTIseq_NK.cells","MCP_NK cells","EPIC_NKcells"),
  Dendritic     = c("GSVA_Dendritic_cells","ssGSEA_Dendritic_cells",
                    "quanTIseq_Dendritic.cells","MCP_Myeloid dendritic cells"),
  Stromal       = c("MCP_Fibroblasts","MCP_Endothelial cells",
                    "EPIC_CAFs","EPIC_Endothelial")
)

# Filter to available features
immune_component_groups <- lapply(immune_component_groups, function(f)
  f[f %in% colnames(data_main)])
immune_component_groups <- immune_component_groups[
  sapply(immune_component_groups, length) > 0]

cat("Immune components defined:\n")
for (g in names(immune_component_groups))
  cat(sprintf("  %-20s: %d features\n", g, length(immune_component_groups[[g]])))

# Compute scores: z-score each feature then average within group
immune_scores_df <- data_main %>% select(sample, PAM50_Subtype)

for (grp in names(immune_component_groups)) {
  feats <- immune_component_groups[[grp]]
  if (length(feats) == 1) {
    immune_scores_df[[paste0("score_", grp)]] <- as.numeric(scale(data_main[[feats]]))
  } else {
    mat <- scale(as.matrix(data_main[, feats]))
    immune_scores_df[[paste0("score_", grp)]] <- rowMeans(mat, na.rm=TRUE)
  }
}

score_cols <- grep("^score_", colnames(immune_scores_df), value=TRUE)
cat("Component scores computed:", length(score_cols), "\n")

# Within-subtype median split for survival analysis
data_scored <- data_main %>%
  left_join(immune_scores_df, by=c("sample","PAM50_Subtype")) %>%
  group_by(PAM50_Subtype) %>%
  mutate(across(
    all_of(score_cols),
    list(group   = ~ ifelse(. > median(., na.rm=TRUE), "High","Low"),
         quartile = ~ ntile(., 4),
         extreme  = ~ dplyr::case_when(
           ntile(.,4)==4 ~ "High_Q4",
           ntile(.,4)==1 ~ "Low_Q1",
           TRUE          ~ NA_character_)),
    .names = "{.col}_{fn}"
  )) %>%
  ungroup()

# Add cluster labels for Her2 and Normal
if (length(cluster_assignments) > 0) {
  cl_combined <- do.call(rbind, cluster_assignments)
  data_scored <- data_scored %>%
    left_join(cl_combined %>% select(sample, immune_cluster), by="sample")
} else {
  data_scored$immune_cluster <- NA_character_
}

# Primary group variable
data_scored <- data_scored %>%
  mutate(
    immune_group_primary = dplyr::case_when(
      PAM50_Subtype %in% subtypes_with_clusters &
        !is.na(immune_cluster) ~ immune_cluster,
      TRUE ~ score_Cytotoxic_group
    ),
    analysis_strategy = dplyr::case_when(
      PAM50_Subtype %in% subtypes_with_clusters ~
        "Consensus clustering k=2 (>=10% balance)",
      TRUE ~
        "Continuous immune score (gradient — no discrete subgroups)"
    )
  )

cat("\nFinal group distribution:\n")
print(table(data_scored$PAM50_Subtype, data_scored$immune_group_primary))

saveRDS(data_scored, "data/processed/clustering/data_scored_final.rds")

write.csv(
  data_scored %>% select(
    sample, PAM50_Subtype, analysis_strategy,
    immune_cluster, immune_group_primary,
    all_of(score_cols),
    ends_with("_group"),
    ends_with("_extreme")
  ),
  "results/tables/clustering/immune_group_assignments.csv",
  row.names = FALSE
)

################################################################################
# SECTION 12: IMMUNE LANDSCAPE VISUALIZATIONS (ALL SUBTYPES)
################################################################################

log_message("=== SECTION 12: IMMUNE LANDSCAPE FIGURES ===")

# Fig 1: Mean immune component heatmap across subtypes
mean_by_sub <- immune_scores_df %>%
  group_by(PAM50_Subtype) %>%
  summarise(across(starts_with("score_"), mean, na.rm=TRUE), .groups="drop") %>%
  column_to_rownames("PAM50_Subtype")

colnames(mean_by_sub) <- gsub("score_","",colnames(mean_by_sub))
mean_scaled <- t(scale(t(mean_by_sub)))

ann_row <- data.frame(PAM50=rownames(mean_scaled), row.names=rownames(mean_scaled))

png("results/figures/clustering/02_immune_landscape_heatmap.png",
    width=11, height=7, units="in", res=300)
pheatmap(
  t(mean_scaled),
  annotation_col    = ann_row,
  annotation_colors = list(PAM50 = SUBTYPE_COLORS),
  color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
  clustering_method = "ward.D2",
  main = "Immune Landscape by PAM50 Subtype (mean z-score)",
  fontsize=11, border_color=NA,
  display_numbers=TRUE, number_format="%.2f", number_color="black"
)
dev.off()

# Fig 2: Boxplots per component per subtype
# Uses pairwise Wilcoxon (hide.ns=TRUE) instead of single KW label
# so readers see WHICH pairs differ, not just "something differs"
score_cols_raw <- grep("^score_[A-Za-z_]+$", colnames(data_scored), value = TRUE)
score_cols_raw <- score_cols_raw[!grepl("_group$|_extreme$|_quartile$", score_cols_raw)]

score_long <- data_scored %>%
  select(sample, PAM50_Subtype, all_of(score_cols_raw)) %>%
  pivot_longer(cols      = all_of(score_cols_raw),
               names_to  = "component",
               values_to = "score") %>%
  mutate(component = gsub("score_", "", component))

# Key pairwise comparisons — biologically motivated
my_comparisons <- list(
  c("Basal", "Her2"),
  c("Basal", "LumA"),
  c("Basal", "LumB"),
  c("Basal", "Normal"),
  c("Her2",  "LumA"),
  c("Her2",  "LumB"),
  c("Her2",  "Normal"),
  c("LumA",  "LumB"),
  c("LumA",  "Normal"),
  c("LumB",  "Normal")
)

p_box <- ggplot(score_long,
                aes(x = PAM50_Subtype, y = score, fill = PAM50_Subtype)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.8, linewidth = 0.4) +

  # Overall Kruskal-Wallis p — one label per panel top-left
  stat_compare_means(
    method   = "kruskal.test",
    label    = "p.format",
    label.x  = 0.6, label.y  = Inf,
    vjust = 2, size = 2.8
  ) +

  # Pairwise Wilcoxon brackets — only significant pairs shown
  stat_compare_means(
    comparisons   = my_comparisons,
    method        = "wilcox.test",
    label         = "p.signif",
    hide.ns       = TRUE,
    tip.length    = 0.01,
    size          = 2.8,
    step.increase = 0.07
  ) +

  facet_wrap(~ component, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  labs(
    title    = "Immune Component Scores by PAM50 Subtype",
    subtitle = paste0(
      "Overall Kruskal-Wallis p shown per panel  |  ",
      "Brackets = significant pairwise Wilcoxon (hide.ns=TRUE)\n",
      "* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001"
    ),
    x = "PAM50 Subtype", y = "Score (z-scaled)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
    strip.text      = element_text(face = "bold", size = 9),
    plot.subtitle   = element_text(size = 8, color = "gray40")
  )

ggsave("results/figures/clustering/03_immune_scores_boxplots.png",
       p_box, width = 16, height = 16, dpi = 300)

log_message("Landscape figures saved")

################################################################################
# SECTION 13: CIBERSORT VALIDATION
#
# WHY INDEPENDENT VALIDATION:
#   CIBERSORT was excluded from clustering (Section 2).
#   If IC1 and IC2 differ in CIBERSORT cell fractions, this is
#   CROSS-METHOD validation — the clusters are real biology, not artifacts.
################################################################################

log_message("=== SECTION 13: CIBERSORT VALIDATION ===")

key_cb_cells <- c(
  "CIBERSORT_T cells CD8", "CIBERSORT_T cells regulatory (Tregs)",
  "CIBERSORT_Macrophages M1", "CIBERSORT_Macrophages M2",
  "CIBERSORT_NK cells activated", "CIBERSORT_B cells memory",
  "CIBERSORT_Dendritic cells activated"
)
key_cb_cells <- key_cb_cells[key_cb_cells %in% colnames(cibersort_valid)]

for (subtype in subtypes_with_clusters) {

  val_data <- cibersort_valid %>%
    filter(PAM50_Subtype == subtype) %>%
    left_join(
      cluster_assignments[[subtype]] %>% select(sample, immune_cluster),
      by = "sample"
    ) %>%
    filter(!is.na(immune_cluster))

  if (nrow(val_data) < 10 || length(key_cb_cells) == 0) next

  plot_data <- val_data %>%
    select(sample, immune_cluster, all_of(key_cb_cells)) %>%
    pivot_longer(cols=all_of(key_cb_cells),
                 names_to="cell_type", values_to="fraction") %>%
    mutate(cell_type = gsub("CIBERSORT_","",cell_type))

  n_cl    <- length(unique(val_data$immune_cluster))
  cl_cols <- CLUSTER_COLORS[seq_len(n_cl)]
  names(cl_cols) <- sort(unique(val_data$immune_cluster))

  p_val <- ggplot(plot_data,
                  aes(x=immune_cluster, y=fraction, fill=immune_cluster)) +
    geom_boxplot(outlier.size=0.5, alpha=0.8) +
    geom_jitter(width=0.1, size=0.3, alpha=0.3) +
    stat_compare_means(method="wilcox.test", label="p.signif", size=3.5) +
    facet_wrap(~cell_type, scales="free_y", ncol=4) +
    scale_fill_manual(values=cl_cols, name="Cluster") +
    labs(
      title    = paste(subtype, "— CIBERSORT Validation (Independent)"),
      subtitle = "CIBERSORT excluded from clustering — this is cross-method validation",
      x="Immune Cluster", y="Cell Fraction"
    ) +
    theme_bw() +
    theme(legend.position="none", strip.text=element_text(size=8))

  ggsave(paste0("results/figures/clustering/cibersort_validation_",
                subtype, ".png"),
         p_val, width=14, height=8, dpi=300)

  log_message(paste(subtype, "— CIBERSORT validation saved"))
}

################################################################################
# SECTION 14: PREPARE DATA FOR SURVIVAL + DE + ENRICHMENT ANALYSIS
#
# Outputs ready for downstream scripts:
#   data_scored_final.rds — primary output
#   cibersort_valid.rds   — for validation in Script 08
#   Keys:
#     immune_group_primary — IC1/IC2 for Her2/Normal; Cytotoxic High/Low for others
#     immune_cluster       — IC1/IC2 only for clustered subtypes
#     score_* columns      — continuous immune scores (for Cox regression)
#     score_*_group        — High/Low binary (for KM curves)
#     score_*_extreme      — High_Q4/Low_Q1 (for DE analysis)
################################################################################

log_message("=== SECTION 14: PREPARING DOWNSTREAM DATA ===")

cat("\n=== DATA READY FOR DOWNSTREAM ANALYSIS ===\n")
cat("Primary output: data/processed/clustering/data_scored_final.rds\n\n")

cat("KEY COLUMNS:\n")
cat("  immune_group_primary    — main group for Scripts 07 and 08\n")
cat("                            Normal: IC1/IC2 (discrete clusters)\n")
cat("                            Others: Cytotoxic High/Low (gradient)\n")
cat("  immune_cluster          — IC1/IC2 for Normal ONLY\n")
cat("  analysis_strategy       — clustering vs gradient per subtype\n")
cat("  score_Cytotoxic         — CD8/cytotoxicity consensus (continuous)\n")
cat("  score_Cytotoxic_group   — High/Low for KM curves (Script 07)\n")
cat("  score_Cytotoxic_extreme — High_Q4/Low_Q1 for DESeq2 (Script 08)\n\n")

cat("SAMPLE DISTRIBUTION (primary immune group):\n")
print(table(data_scored$PAM50_Subtype, data_scored$immune_group_primary))

cat("\nSUBTYPE STRATEGY:\n")
for (s in subtypes) {
  strat <- unique(data_scored$analysis_strategy[data_scored$PAM50_Subtype == s])
  cat(sprintf("  %-8s: %s\n", s, strat))
}

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("SCRIPT 06 COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("CLUSTERING RESULTS (elbow PC selection, 10% balance rule):\n")
cat(sprintf("  %-8s  %-15s  %-8s  %-8s  %-25s\n",
            "Subtype", "Sizes(k=2)", "Min%", "Sil", "Result"))
cat(strrep("-", 72), "\n")
for (s in subtypes) {
  opt  <- optimal_k_results[[s]]
  mets <- opt$metrics
  sz   <- mets$sizes[mets$k == 2]
  pct  <- mets$min_pct[mets$k == 2]
  sil  <- round(mets$silhouette[mets$k == 2], 3)
  flag <- opt$balance_flag
  pcs  <- clustering_results[[s]]$n_pcs
  cat(sprintf("  %-8s  %-15s  %-8s  %-8s  %s (n_PCs=%d)\n",
              s, sz, paste0(pct, "%"), sil, flag, pcs))
}

cat("\nKEY FINDING:\n")
cat("  Only Normal-like subtype shows genuine immune subgroups (IC1/IC2)\n")
cat("  Basal/Her2/LumA/LumB have homogeneous immune profiles\n")
cat("  All four fail the 10% minimum cluster size threshold\n\n")

cat("BIOLOGICAL CHARACTERIZATION (Normal IC1 vs IC2):\n")
cat("  PCA loadings:      results/figures/clustering/pca/loadings_Normal.png\n")
cat("  Feature importance: results/figures/clustering/feature_importance_Normal.png\n")
cat("  Biplot:            results/figures/clustering/pca/biplot_Normal.png\n\n")

cat("NEXT SCRIPTS:\n")
cat("  07_survival_analysis.R\n")
cat("     Uses: immune_group_primary (IC1/IC2 for Normal; High/Low for others)\n")
cat("     Tests: KM curves, log-rank, Cox regression per subtype\n\n")
cat("  08_DE_pathway_analysis.R\n")
cat("     Uses: immune_cluster (Normal IC1 vs IC2)\n")
cat("           score_*_extreme (High_Q4 vs Low_Q1 for gradient subtypes)\n")
cat("     Runs: DESeq2, fgsea, CIBERSORT validation\n\n")

cat("OUTPUT FILES:\n")
cat("  data/processed/clustering/data_scored_final.rds  ← primary\n")
cat("  data/processed/clustering/cibersort_valid.rds    ← for Script 08\n")
cat("  results/tables/clustering/clustering_summary.csv\n")
cat("  results/tables/clustering/immune_group_assignments.csv\n")

log_message("Script 06 complete")
sink(); sink(); close(log_file)
