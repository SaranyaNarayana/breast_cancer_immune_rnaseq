############################################################
# Normalization, filtering, and quality control 
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - 
#   - 
#   - 


############################################################
setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")
############################################################

#-----------------------------------
# Logging
#-----------------------------------

log_file <- file(
  paste0("logs/05_immune_deconvolution_",
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
    library(GSVA)
    library(immunedeconv)
    #library(ESTIMATE)
    library(GSEABase)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(pheatmap)
    library(msigdbr)
})

#create output directory
dir.create("results/figures/immune", recursive=TRUE, showWarnings = FALSE)
dir.create("data/processed/immune", recursive=TRUE, showWarnings = FALSE)


## -----------------------------
## 1. Load Data
## -----------------------------
expr <- readRDS("data/processed/expression_vst_normalized.rds") #  DESeq2 VST as primary normalized count data 
cat("Expression data loaded. Dimensions:", dim(expr), "\n")



clinical <- readRDS("data/processed/clinical_filtered.rds")
cat("Clinical data loaded. Dimensions:", dim(clinical), "\n")


## -----------------------------
## 2. Define immune gene signatures
## -----------------------------

#1. Hallmark immune signatures from MSigDB

library(msigdbr)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
    filter(grepl("IMMUNE|INFLAMMATORY|INTERFERON|TNF|\\bIL[0-9_]", gs_name))