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
    total_counts=colSums(counts), #totoal counts per sample
    n_genes_detected=colSums(counts>0) #unique genes per sample that have at least one read.
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


#Set thresholds for filtering
qc_thresholds <- list(
  # Total counts
  min_total_counts = 15e6, # Minimum total counts (15 million)
  max_total_counts = Inf, # No upper limit for total counts
  
  # Genes detected  
  min_genes_detected = 20000, 
  max_genes_detected =45000  
)

# Apply filters
severe_outliers <- qc_metrics_sample %>%
  filter(
    total_counts < qc_thresholds$min_total_counts |
    total_counts > qc_thresholds$max_total_counts |
    n_genes_detected < qc_thresholds$min_genes_detected |
    n_genes_detected > qc_thresholds$max_genes_detected
  )

cat("Number of severe outliers identified:", nrow(severe_outliers), "\n")#1

samples_to_keep <- setdiff(colnames(counts), severe_outliers$sample)

# Filter count matrices

counts_clean <- counts[, samples_to_keep]
tpm_clean <- tpm[, samples_to_keep]

# Filter clinical data

clinical_clean <- clinical %>%
  filter(barcode %in% samples_to_keep)

cat("Samples after QC:", ncol(counts_clean), "\n") #1013
cat("Retention rate:",round(ncol(counts_clean)/ncol(counts) * 100, 1), "%\n")


## -----------------------------
## 3. Quality control-Genes
## -----------------------------

# Calculate gene-level metrics based on counts
qc_metrics_genes <- data.frame(
  gene=rownames(counts_clean),
  mean_expression=rowMeans(counts_clean), # mean expression across all samples for each gene
  n_samples_expressed=rowSums(counts_clean > 0),# number of samples where the gene is expressed (count > 0)
  pct_samples_expressed=rowSums(counts_clean > 0) / ncol(counts_clean) * 100
)

cat("Summary of qc_metrics_genes:\n")
summary(qc_metrics_genes)


#Caculate gene-level metrics based on TPM
Qc_metircs_tpm <- data.frame(
  gene=rownames(tpm_clean),
  mean_tpm=rowMeans(tpm_clean), # mean TPM across all samples for each gene
  n_samples_expressed=rowSums(tpm_clean > 1),# number of samples where the gene is expressed (TPM > 1)
  pct_samples_expressed=rowSums(tpm_clean > 1) / ncol(tpm_clean) * 100
)
cat("Summary of Qc_metircs_tpm:\n")
summary(Qc_metircs_tpm)

p3<- ggplot(qc_metrics_genes, aes(x=pct_samples_expressed))+
    geom_histogram(bins=50, fill="lightgreen", alpha=0.7)+
    labs(title="Percentage of samples with gene expression",
    x="Percentage of samples gene expressed (TPM>1)",
    y="Number of genes")+
    theme_bw()

ggsave("results/figures/preprocessing_plots/03_Percentage_of_samples_with_gene_expression_TPM.png", p3, width = 8, height = 6)


#filter the low expressed genes based on TPM: keep genes with TPM >1 in at least 10% of the samples 
min_samples_expressed <- ceiling(0.10 * ncol(tpm_clean)) #10% of samples; 102
cat("Minimum number of samples a gene must be expressed in to be retained:", min_samples_expressed, "\n")#102

genes_to_keep <- rowSums(tpm_clean > 1) >= min_samples_expressed #genes that have TPM >1 in at least 10% of samples

cat("Number of genes retained after filtering:", sum(genes_to_keep), "\n") # 23059
cat("Retention rate:", round(sum(genes_to_keep) / nrow(counts_clean) * 100, 1), "%\n")#38%


#Extract the counts and TPM for the retained genes, and the corresponding clinical data
counts_filtered <- counts_clean[genes_to_keep, ]
tpm_filtered <- tpm_clean[genes_to_keep, ]
clinical_filtered <- clinical_clean

cat("Dimensions of filtered counts matrix:", dim(counts_filtered), "\n") #23059 1013
cat("Dimensions of filtered TPM matrix:", dim(tpm_filtered), "\n") #23059 1013
cat("Dimensions of filtered clinical data:", dim(clinical_filtered), "\n") #1013 19

saveRDS(counts_filtered, "data/processed/counts_filtered.rds")
saveRDS(tpm_filtered, "data/processed/tpm_filtered.rds")
saveRDS(clinical_filtered, "data/processed/clinical_filtered.rds")

# counts_filtered <- readRDS("data/processed/counts_filtered.rds")
# tpm_filtered <- readRDS("data/processed/tpm_filtered.rds")
# clinical_filtered <- readRDS("data/processed/clinical_filtered.rds")

## -----------------------------
## 4. Normalization for downstream analaysis
## -----------------------------


# Method 1: TMM normalization (edgeR)
dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge, method = "TMM")
tmm_norm <- cpm(dge, log = TRUE, prior.count = 1)

# Method 2: DESeq2 normalization (size factors)
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = clinical_filtered,
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
deseq2_norm <- vst(dds, blind = TRUE)  # Variance stabilizing transformation

# Method 3: Log2(TPM + 1)
log2_tpm <- log2(tpm_filtered + 1)

# Save normalized data
saveRDS(tmm_norm, "data/processed/expression_tmm_normalized.rds")
saveRDS(deseq2_norm, "data/processed/expression_vst_normalized.rds")
saveRDS(log2_tpm, "data/processed/expression_log2tpm.rds")

# Use DESeq2 VST as primary normalized data for downstream analysis
expr_norm <- assay(deseq2_norm)