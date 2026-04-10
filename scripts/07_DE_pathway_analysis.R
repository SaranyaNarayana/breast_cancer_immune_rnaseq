################################################################################
# 08_DE_pathway_analysis.R
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
#
# PURPOSE:
#   Identify differentially expressed genes (DEGs) between immune subgroups
#   and determine which biological pathways are enriched in each group.
#
# RESEARCH QUESTION:
#   What genes and biological processes distinguish immune subgroups
#   within each PAM50 subtype?
#
# ANALYSIS STRATEGY PER SUBTYPE:
#   Normal:              IC1 vs IC2 (discrete clusters from Script 06)
#   Basal/Her2/LumA/LumB: Cytotoxic High_Q4 vs Low_Q1 (continuous gradient)
#
# WHY Q4 vs Q1 (NOT HIGH vs LOW):
#   Median split (High/Low) includes all samples — the middle is noisy.
#   Q4 (top 25%) vs Q1 (bottom 25%) compares the EXTREME ends of the gradient.
#   This maximises signal-to-noise and gives cleaner DEGs.
#   Trade-off: smaller n, but for DE signal quality > sample size.
#
# PIPELINE:
#   Script 06 data_scored_final.rds → THIS SCRIPT → Script 07
#
# OUTPUTS:
#   results/tables/DE/         — DESeq2 results per subtype
#   results/tables/pathway/    — fgsea results per subtype
#   results/figures/DE/        — volcano plots, MA plots, heatmaps
#   results/figures/pathway/   — pathway bubble plots
#   data/processed/DE/         — R objects for downstream use
################################################################################

setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")

log_file <- file(
  paste0("logs/07_DE_pathway_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"),
  open = "wt"
)
sink(log_file, type = "output")
sink(log_file, type = "message")

log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(cowplot)
  library(msigdbr)
})


# Output directories
dir.create("results/tables/DE",       recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pathway",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/DE",      recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/pathway", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/DE",       recursive = TRUE, showWarnings = FALSE)

SUBTYPE_COLORS <- c(
  Basal="#E41A1C", Her2="#FF7F00",
  LumA="#4DAF4A", LumB="#377EB8", Normal="#984EA3"
)
CLUSTER_COLORS <- c(IC1="#E41A1C", IC2="#377EB8")

################################################################################
# SECTION 1: LOAD DATA
#
# WHY TWO SOURCES:
#   data_scored_final.rds  — group assignments from Script 06
#   raw_counts.rds         — DESeq2 REQUIRES raw integer counts, not normalised
#                            values. TPM/FPKM/log-normalised data gives wrong results.
################################################################################

log_message("=== SECTION 1: LOADING DATA ===")

data_scored <- readRDS("data/processed/clustering/data_scored_final.rds")
cat("data_scored dimensions:", dim(data_scored), "\n")#1013 192 
colnames(data_scored)
library(dplyr)

data_scored %>%
  select(sample, GSVA_T_cells_CD8, GSVA_Cytotoxicity,
         ssGSEA_T_cells_CD8, ssGSEA_Cytotoxicity,
         quanTIseq_T.cells.CD8, `MCP_CD8 T cells`, EPIC_CD8_Tcells,
         score_Cytotoxic, score_Cytotoxic_group, 
         score_Cytotoxic_quartile, score_Cytotoxic_extreme) %>%
  head() %>%
  as.data.frame()



required_cols <- c("sample","PAM50_Subtype","immune_cluster",
                   "immune_group_primary","score_Cytotoxic_extreme")
head(data_scored [, required_cols])

missing <- setdiff(required_cols, colnames(data_scored))
if (length(missing) > 0)
  stop("Missing columns: ", paste(missing, collapse=", "),
       "\nRe-run Script 06 first.")

##--------------------------------------------------------------
# Raw count matrix- genes x samples, raw integer counts
##-------------------------------------------------------------
count_matrix <- readRDS("data/processed/counts_filtered.rds")
# Alternative: count_matrix <- read.csv("data/processed/expression/raw_counts.csv", row.names=1)

cat("Count matrix dimensions:", dim(count_matrix), "\n") # 23059 1013 
cat("First gene names:", head(rownames(count_matrix), 5), "\n")
cat("First sample names:", head(colnames(count_matrix), 5), "\n")
head(count_matrix[, 1:5])



# Ensure integer type (DESeq2 requires this)
if (!is.integer(count_matrix[1,1])) {
  mode(count_matrix) <- "integer"
  cat("Converted to integer\n")
}


################################################################################
# SECTION 1b: ENSEMBL ID → GENE SYMBOL MAPPING
log_message("=== SECTION 1b: ENSEMBL TO SYMBOL MAPPING ===")

id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rownames(count_matrix),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

cat("Raw ID mapping dimensions:", dim(id_map), "\n")
cat("Unique Ensembl IDs in map:", length(unique(id_map$ENSEMBL)), "\n")
cat("Unmapped (NA symbol):", sum(is.na(id_map$SYMBOL)), "\n")

# Step 2: Remove unmapped and duplicate symbols
id_map_clean <- id_map[
  !is.na(id_map$SYMBOL) & !duplicated(id_map$SYMBOL),
]

cat("Filtered ID mapping dimensions:", dim(id_map_clean), "\n")
cat("Unique symbols retained:", length(unique(id_map_clean$SYMBOL)), "\n")


# Step 3: Subset count matrix to only the mapped Ensembl IDs
# Reorder rows to exactly match id_map_clean (required for correct symbol assignment)
count_matrix_symbol <- count_matrix[id_map_clean$ENSEMBL, ]

cat("Symbol-mapped count matrix dimensions:", dim(count_matrix_symbol), "\n")

# Step 4: Replace Ensembl IDs with gene symbols as row names
rownames(count_matrix_symbol) <- id_map_clean$SYMBOL

# Step 5: Verify no duplicate rownames remain
n_dup <- sum(duplicated(rownames(count_matrix_symbol)))
cat("Duplicate gene symbols remaining:", n_dup, "\n")
if (n_dup > 0) stop("Duplicate rownames found — check id_map filtering step.")

cat("First 5 gene symbols:\n")
print(head(rownames(count_matrix_symbol), 5))

# Step 6: Confirm integer type is preserved
if (!is.integer(count_matrix_symbol[1,1])) {
  mode(count_matrix_symbol) <- "integer"
}

# Step 7: Save mapped matrix for reuse
# Loading this avoids re-running the mapping on subsequent runs
saveRDS(count_matrix_symbol,
        "data/processed/DE/count_matrix_symbol.rds")


# Check sample overlap between count matrix and data_scored
overlap <- intersect(colnames(count_matrix_symbol), data_scored$sample)
cat("Sample overlap (count matrix vs data_scored):", length(overlap), "\n")
cat("Samples in counts NOT in data_scored:",
    length(setdiff(colnames(count_matrix_symbol), data_scored$sample)), "\n")
cat("Samples in data_scored NOT in counts:",
    length(setdiff(data_scored$sample, colnames(count_matrix_symbol))), "\n")

if (length(overlap) < 100)
  warning("Low sample overlap — check sample name format in both objects.")

# Replace count_matrix with symbol version for all downstream use
count_matrix <- count_matrix_symbol
cat("\ncount_matrix now uses gene symbols as rownames.\n")
cat("Final dimensions:", dim(count_matrix), "\n\n") #18568 1013



subtypes <- c("Basal","Her2","LumA","LumB","Normal")




################################################################################
# SECTION 2: DEFINE COMPARISON GROUPS
#
# Normal   : IC1 (immune-cold) vs IC2 (immune-active) — discrete clusters
# All other: Cytotoxic Low_Q1 (ref) vs High_Q4 (comp) — gradient extremes
#
# WHY CYTOTOXIC SCORE AS RANKING AXIS:
#   Most biologically meaningful and consistent immune axis across subtypes.
#   Integrates CD8 signal from 5 independent methods.
#   Change to score_Inflammatory_extreme etc if preferred.
################################################################################

log_message("=== SECTION 2: DEFINING COMPARISON GROUPS ===")

get_comparison_groups <- function(data_scored, subtype) {

  sub_data <- data_scored %>% filter(PAM50_Subtype == subtype)

  if (subtype == "Normal") {
    groups <- sub_data %>%
      filter(!is.na(immune_cluster)) %>%
      dplyr::select(sample, group = immune_cluster)
    ref   <- "IC1"
    comp  <- "IC2"
    label <- "IC1_vs_IC2"
  } else {
    groups <- sub_data %>%
      filter(score_Cytotoxic_extreme %in% c("High_Q4","Low_Q1")) %>%
      dplyr::select(sample, group = score_Cytotoxic_extreme)
    ref   <- "Low_Q1"
    comp  <- "High_Q4"
    label <- "High_Q4_vs_Low_Q1"
  }

  cat(sprintf("  %s | %s | n=%d\n", subtype, label, nrow(groups)))
  print(table(groups$group))
  return(list(groups=groups, ref=ref, comp=comp, label=label))
}

comparison_list <- list()
for (subtype in subtypes) {
  cat("\n", subtype, ":\n", sep="")
  comparison_list[[subtype]] <- get_comparison_groups(data_scored, subtype)
}



################################################################################
# SECTION 3: DESeq2 DIFFERENTIAL EXPRESSION
#
# DESeq2 models RNA-seq counts with a negative binomial distribution.
# This accounts for:
#   (a) Non-negative integer nature of counts
#   (b) Mean-variance relationship (overdispersion)
#   (c) Library size differences between samples (size factor normalisation)
#
# PIPELINE STEPS:
#   1. estimateSizeFactors  — library size normalisation
#   2. estimateDispersions  — gene-wise dispersion (biological variability)
#   3. nbinomWaldTest       — Wald test per gene for group difference
#   4. lfcShrink(apeglm)    — shrink noisy LFCs for low-count genes
#
# THRESHOLDS:
#   |log2FC| > 1   — at least 2-fold change
#   padj < 0.05    — BH-corrected FDR
#   baseMean > 10  — exclude very lowly expressed genes
################################################################################

log_message("=== SECTION 3: DESeq2 DIFFERENTIAL EXPRESSION ===")

run_deseq2 <- function(count_matrix, groups_df, ref, comp,
                        subtype, label,
                        lfc_threshold  = 1.0,
                        padj_threshold = 0.05,
                        basemean_min   = 10) {

  cat("\n--- DESeq2:", subtype, "|", label, "---\n")

  common_samples <- intersect(colnames(count_matrix), groups_df$sample)
  cat("Common samples:", length(common_samples), "\n")

  if (length(common_samples) < 6) {
    cat("WARNING: Too few samples — skipping\n")
    return(NULL)
  }

  counts_sub <- count_matrix[, common_samples, drop=FALSE]
  meta_sub   <- groups_df %>%
    filter(sample %in% common_samples) %>%
    arrange(match(sample, common_samples)) %>%
    mutate(group = factor(group, levels=c(ref, comp)))

  cat("Group sizes:\n"); print(table(meta_sub$group))

  # Build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = meta_sub,
    design    = ~ group
  )

  # Pre-filter: at least 5 counts in at least 10% of samples
  keep <- rowSums(counts(dds) >= 5) >= max(3, floor(ncol(dds)*0.10))
  dds  <- dds[keep, ]
  cat("Genes after pre-filter:", nrow(dds), "\n")

  # Run DESeq2
  dds <- DESeq(dds, quiet=TRUE)

  # LFC shrinkage — apeglm is the recommended posterior estimation method
  # It shrinks unreliable LFCs toward 0 without affecting significant ones
  coef_name <- paste0("group_", comp, "_vs_", ref)
  cat("Shrinking LFC:", coef_name, "\n")

  res_shrunk <- lfcShrink(dds, coef=coef_name, type="apeglm", quiet=TRUE)

  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene") %>%
    mutate(
      subtype    = subtype,
      comparison = label,
      direction  = dplyr::case_when(
        log2FoldChange >  lfc_threshold & padj < padj_threshold ~ "UP",
        log2FoldChange < -lfc_threshold & padj < padj_threshold ~ "DOWN",
        TRUE                                                     ~ "NS"
      ),
      significant = direction != "NS"
    ) %>%
    filter(!is.na(padj), baseMean >= basemean_min) %>%
    arrange(padj)

  n_up   <- sum(res_df$direction=="UP")
  n_down <- sum(res_df$direction=="DOWN")
  cat(sprintf("  DEGs: %d UP | %d DOWN\n", n_up, n_down))

  write.csv(res_df,
    paste0("results/tables/DE/DESeq2_", subtype, "_", label, ".csv"),
    row.names=FALSE)

  return(list(res=res_df, dds=dds, n_up=n_up, n_down=n_down,
               subtype=subtype, label=label))
}

de_results <- list()
for (subtype in subtypes) {
  cmp <- comparison_list[[subtype]]
  de_results[[subtype]] <- run_deseq2(
    count_matrix, cmp$groups, cmp$ref, cmp$comp,
    subtype, cmp$label
  )
}

saveRDS(de_results, "data/processed/DE/de_results_all.rds")

##------------------------------------------------------------------##
#Quality Check
# Run this to verify the direction is correct
# UP should contain immune genes, DOWN should contain tumour/stromal genes
##------------------------------------------------------------------##
for (subtype in subtypes) {
  res <- de_results[[subtype]]$res

  cat("\n===", subtype, "===\n")

  # Top 10 UP genes
  cat("Top 10 UP (higher in immune-active):\n")
  print(res %>% filter(direction=="UP") %>%
          head(10) %>% dplyr::select(gene, log2FoldChange, padj))

  # Top 10 DOWN genes
  cat("Top 10 DOWN (higher in immune-cold):\n")
  print(res %>% filter(direction=="DOWN") %>%
          head(10) %>% dplyr::select(gene, log2FoldChange, padj))
}







################################################################################
# SECTION 4: VOLCANO PLOTS
#
# X = log2FC (effect size)   Y = -log10(padj) (significance)
# Top-right = UP + significant    Top-left = DOWN + significant
# Top genes labelled with ggrepel (no overlap)
################################################################################

log_message("=== SECTION 4: VOLCANO PLOTS ===")

plot_volcano <- function(res_df, subtype, label,
                          lfc_cut=1.0, padj_cut=0.05, n_label=15) {

  if (is.null(res_df) || nrow(res_df)==0) return(NULL)

  plot_df <- res_df %>%
    mutate(log10_padj = -log10(padj + 1e-300))

  top_genes <- plot_df %>% filter(significant) %>%
    arrange(padj) %>% head(n_label)

  n_up   <- sum(plot_df$direction=="UP")
  n_down <- sum(plot_df$direction=="DOWN")

  p <- ggplot(plot_df, aes(x=log2FoldChange, y=log10_padj,
                            color=direction)) +
    geom_point(alpha=0.5, size=1.2) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut),
               linetype="dashed", color="gray40", linewidth=0.6) +
    geom_hline(yintercept=-log10(padj_cut),
               linetype="dashed", color="gray40", linewidth=0.6) +
    geom_text_repel(data=top_genes, aes(label=gene),
                    size=2.8, max.overlaps=20, fontface="italic") +
    annotate("text", x=max(plot_df$log2FoldChange)*0.7,
             y=max(plot_df$log10_padj)*0.95,
             label=paste0(n_up," UP"), color="#E41A1C", fontface="bold", size=4) +
    annotate("text", x=min(plot_df$log2FoldChange)*0.7,
             y=max(plot_df$log10_padj)*0.95,
             label=paste0(n_down," DOWN"), color="#377EB8", fontface="bold", size=4) +
    scale_color_manual(values=c("UP"="#E41A1C","DOWN"="#377EB8","NS"="gray75"),
                       name=NULL) +
    labs(
      title    = paste(subtype, "— Volcano Plot"),
      subtitle = paste0(label, "\n|log2FC|>", lfc_cut, "  padj<", padj_cut),
      x="log2 Fold Change", y="-log10(adjusted p-value)"
    ) + theme_bw() +
    theme(legend.position="bottom",
          plot.subtitle=element_text(size=8, color="gray40"))

  ggsave(paste0("results/figures/DE/volcano_", subtype, "_", label, ".png"),
         p, width=9, height=7, dpi=300)
  return(p)
}

for (subtype in subtypes) {
  if (!is.null(de_results[[subtype]]))
    plot_volcano(de_results[[subtype]]$res, subtype,
                 comparison_list[[subtype]]$label)
}




################################################################################
# SECTION 5: MA PLOTS
#
# X = log10(baseMean) — average expression
# Y = log2FC — fold change
# Checks: is DE biased toward low/high expressed genes?
# LFC shrinkage should narrow the "funnel" at low expression
################################################################################

log_message("=== SECTION 5: MA PLOTS ===")

plot_ma <- function(res_df, subtype, label) {
  if (is.null(res_df) || nrow(res_df)==0) return(NULL)

  plot_df <- res_df %>% mutate(log_bm = log10(baseMean+1))

  ggplot(plot_df, aes(x=log_bm, y=log2FoldChange, color=direction)) +
    geom_point(alpha=0.4, size=1) +
    geom_hline(yintercept=0, linewidth=0.8) +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="gray50") +
    scale_color_manual(values=c("UP"="#E41A1C","DOWN"="#377EB8","NS"="gray75"),
                       name=NULL) +
    labs(title=paste(subtype,"— MA Plot (apeglm LFC shrinkage)"),
         subtitle=label, x="log10(baseMean)", y="log2FC") +
    theme_bw()

  ggsave(paste0("results/figures/DE/MA_", subtype, "_", label, ".png"),
         last_plot(), width=8, height=6, dpi=300)
}

for (subtype in subtypes)
  if (!is.null(de_results[[subtype]]))
    plot_ma(de_results[[subtype]]$res, subtype,
            comparison_list[[subtype]]$label)




