############################################################
# Normalization, filtering, and quality control 
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - 
#   - 


############################################################
setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")
############################################################

#-----------------------------------
# Logging
#-----------------------------------

log_file <- file(
  paste0("logs/04_preprocessing_",
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
    library(DESeq2)
    library(edgeR)
    library(limma)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
})

#install.packages("pheatmap")


# Create output directories
dir.create("results/figures/preprocessing_plots", recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 1. Load Data
## -----------------------------
counts <-readRDS("data/processed/counts_match.rds")
cat("counts", dim(counts), "\n") #60660 1014

tpm <- readRDS("data/processed/tpm_unstrand_match.rds")
cat("tpm", dim(tpm), "\n") #60660 1014

clinical <- readRDS("data/processed/clinical_filtered_clean.rds")
cat("clinical data", dim(clinical), "\n") #1014 19 

 

## -----------------------------
## 2. Quality control-Samples
## -----------------------------

#Total count distribution

qc_metrics_sample <- data.frame(
    sample=colnames(counts),
    total_counts=colSums(counts),
    n_genes_detected=colSums(counts>0) #counts how many unique genes in a sample that have at least one read.
)
cat("qc_metrics_sample", dim(qc_metrics_sample), "\n")
head(qc_metrics_sample )

cat("\nSummary of qc_metrics_sample:\n")
summary(qc_metrics_sample)

# plot 1: Total count distribution
p1 <- ggplot (qc_metrics_sample, aes(x=total_counts/1e6))+
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(title = "Total Read Counts per Sample",
       x = "Total Counts (millions)",
       y = "Number of Samples") +
  theme_bw()
ggsave("results/figures/preprocessing_plots/01_total_counts_distribution.png", p1, width = 8, height = 6)

# Plot 2: genes detected

p2 <- ggplot(qc_metrics_sample, aes(x=n_genes_detected))+
    geom_histogram(bins=50, fill="coral", alpha=0.7)+
    labs(title="Genes detected per sample",
    x="Number of genes detected",
    y="Number of sample")+
    theme_bw()
ggsave("results/figures/preprocessing_plots/02_genes_detected_distribution.png", p2, width = 8, height = 6)


# Look at quantiles
quantile(qc_metrics_sample$total_counts / 1e6, probs = c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99))
quantile(qc_metrics_sample$n_genes_detected, probs = c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99))


#Set thresholds for filtering based on quantiles
qc_thresholds <- list(
  # Total counts
  min_total_counts = quantile(qc_metrics_sample$total_counts, 0.01),
  max_total_counts = quantile(qc_metrics_sample$total_counts, 0.99),    
  
  # Genes detected  
  min_genes_detected = quantile(qc_metrics_sample$n_genes_detected, 0.01), 
  max_genes_detected = quantile(qc_metrics_sample$n_genes_detected, 0.99)  
)

# Apply filters
severe_outliers <- qc_metrics_sample %>%
  filter(
    total_counts < qc_thresholds$min_total_counts |
    total_counts > qc_thresholds$max_total_counts |
    n_genes_detected < qc_thresholds$min_genes_detected |
    n_genes_detected > qc_thresholds$max_genes_detected
  )

cat("Number of severe outliers identified:", nrow(severe_outliers), "\n")#38

samples_to_keep <- setdiff(colnames(counts), severe_outliers$sample)

# Filter count matrices

counts_clean <- counts[, samples_to_keep]
tpm_clean <- tpm[, samples_to_keep]

# Filter clinical data

clinical_clean <- clinical %>%
  filter(barcode %in% samples_to_keep)

cat("Samples after QC:", ncol(counts_clean), "\n")
cat("Retention rate:",
round(ncol(counts_clean)/ncol(counts) * 100, 1), "%\n")**