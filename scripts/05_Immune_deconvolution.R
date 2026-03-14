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
    library(dplyr)
    library(limma)
    library(biomaRt)
})

#create output directory
dir.create("results/figures/immune", recursive=TRUE, showWarnings = FALSE)
dir.create("data/processed/immune", recursive=TRUE, showWarnings = FALSE)
dir.create("data/immune", recursive=TRUE, showWarnings = FALSE)


## -----------------------------
## 1. Load Data
## -----------------------------
expr <- readRDS("data/processed/expression_vst_normalized.rds") #  DESeq2 VST as primary normalized count data 
cat("Expression data loaded. Dimensions:", dim(expr), "\n")

tpm <- readRDS("data/processed/tpm_filtered.rds")
cat("TPM dimensions:", dim(tpm), "\n")

clinical <- readRDS("data/processed/clinical_filtered.rds")
cat("Clinical data loaded. Dimensions:", dim(clinical), "\n")


## -----------------------------
## 2. Define immune gene signatures
## -----------------------------

# Hallmark immune signatures from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    filter(grepl("IMMUNE|INFLAMMATORY|INTERFERON|TNF|\\bIL[0-9_]", gs_name))
cat("Hallmark immune sets dimensions:", dim(hallmark_sets), "\n")

saveRDS(hallmark_sets, "data/immune/hallmark_sets.rds")

write.csv(hallmark_sets, "data/immune/hallmark_sets.csv", row.names = FALSE)


#Convert to list format for GSVA

hallmark_list <- hallmark_sets %>% 
    split(.$gs_name)%>%
    lapply(function(x) x$gene_symbol)

cat("Hallmark immune signatures:", length(hallmark_list), "\n")#4


# Custom immune signatures for literature
immune_signatures <- list(
  
  # T cell markers
  T_cells_CD8 = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "EOMES", "TBX21"),
  T_cells_CD4 = c("CD4", "CD40LG", "ICOS", "IL2", "IL21"),
  T_reg = c("FOXP3", "IL2RA", "IKZF2", "CTLA4", "TIGIT"),
  T_exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "BTLA", "CD244"),
  
  # B cells and plasma cells
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "BLK", "BANK1"),
  Plasma_cells = c("SLAMF7", "TNFRSF17", "SDC1", "TNFRSF13B"),
  
  # Myeloid cells
  Macrophages_M1 = c("NOS2", "IRF5", "CD80", "CD86", "IL1B", "IL12A", "TNF"),
  Macrophages_M2 = c("CD163", "MRC1", "MSR1", "IL10", "ARG1", "TGFB1"),
  Dendritic_cells = c("CD1C", "CLEC9A", "XCR1", "BATF3", "IRF8", "FLT3"),
  
  # NK cells
  NK_cells = c("KLRB1", "KLRK1", "NCR1", "NKG7", "GNLY", "FCGR3A"),
  
  # Immune checkpoints
  Checkpoints = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", 
                  "TIGIT", "BTLA", "VSIR"),
  
  # Cytokines and chemokines
  Chemokines = c("CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4", "CCL5"),
  Cytokines_inflammatory = c("IFNG", "IL2", "IL12A", "IL12B", "TNF", "IL1B"),
  
  # Functional signatures
  Cytotoxicity = c("PRF1", "GZMA", "GZMB", "GZMH", "GNLY", "NKG7", "KLRK1"),
  Antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", 
                           "B2M", "TAP1", "TAP2"),
  IFN_gamma_response = c("IFNG", "IRF1", "STAT1", "IDO1", "CXCL9", "CXCL10", 
                         "CXCL11", "GBP1"),
  IFN_alpha_response = c("IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "ISG15"),
  
  # Immunosuppressive
  Immunosuppression = c("TGFB1", "IL10", "VEGFA", "IDO1", "ARG1", "CD274")
)

# Combine all signatures list
all_signatures <- c(hallmark_list, immune_signatures)

cat("\nTotal immune signatures:", length(all_signatures), "\n")#22

saveRDS(all_signatures, "data/immune/all_signatures.rds")



## ----------------------------------------------------------------------------------------------
## 3. Convert Ensembl IDs in the  expr_matrix to gene symbols using AnnotationDbi and org.Hs.eg.db
## ----------------------------------------------------------------------------------------------
# Convert expression to matrix 
head(expr)[,1:5] # has Ensembl IDs

head(all_signatures[[1]]) #has gene symbols, need to convert to Ensembl IDs

#Convert Ensembl IDs in the  expr_matrix to gene symbols using AnnotationDbi and org.Hs.eg.db
install.packages(c("AnnotationDbi", "org.Hs.eg.db"))
library(AnnotationDbi)
library(org.Hs.eg.db)  # install via BiocManager::install("org.Hs.eg.db")

# Strip version suffix from rownames
rownames(expr) <- sub("\\..*", "", rownames(expr_matrix))
head(rownames(expr))


#Extract gene symbols of Ensembl IDs of the expression matrix using AnnotationDbi 
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rownames(expr),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)
head(id_map)
cat("ID mapping dimensions:", dim(id_map), "\n") #23267 2 

# Remove unmapped or duplicate symbols (Keep only first symbol per Ensembl ID)
id_map_1 <- id_map[!is.na(id_map$SYMBOL) & !duplicated(id_map$SYMBOL), ]
cat("Filtered ID mapping dimensions:", dim(id_map_1), "\n") # 18568 2 


# Update matrix rownames
expr_filtered_symbol <- expr[id_map_1$ENSEMBL, ] #Reorder to match id_map_1 exactly
cat("Filtered expression filtered symbol dimensions:", dim(expr_filtered_symbol), "\n") # 18568 1013

rownames(expr_filtered_symbol) <- id_map_1$SYMBOL #Assign gene symbols
head(rownames(expr_filtered_symbol))

saveRDS(expr_filtered_symbol, "data/processed/immune/expression_filtered_symbol.rds")

# Check how many of your lost genes overlap with your gene sets
all_set_genes <- unique(unlist(all_signatures))
length(all_set_genes) # 574


# Genes in sets that are still in your matrix
covered <- intersect(rownames(expr_filtered_symbol), all_set_genes)

cat("Genes required by sets:", length(all_set_genes), "\n")# 574
cat("Covered after mapping:", length(covered), "\n") # 552
cat("Coverage (%): ", round(length(covered)/length(all_set_genes)*100, 1), "%\n") # 96.2%


# Check size of each gene set before gsvaParam filtering to decide the thresholds for minSize and maxSize
set_sizes <- sapply(all_signatures, length)
print(sort(set_sizes))

# Check effective size (overlap with expression matrix)
effective_sizes <- sapply(all_signatures, function(gs) {
  length(intersect(gs, rownames(expr_filtered_symbol)))
})

# Summary
data.frame(
  gene_set       = names(effective_sizes),
  original_size  = sapply(all_signatures, length),
  effective_size = effective_sizes
) |> print()



## -----------------------------
## 4. GENE SET ENRICHMENT ANALYSIS
## -----------------------------

## 4a. Gene set variation analysis (Scores -1 to 1 relative to other samples);
##     Which samples have relatively higher or lower pathway activity compared to others?

cat("\nRunning GSVA...\n")

expr_matrix <- as.matrix(expr_filtered_symbol)

# Build the parameter object for GSVA 
gsva_scores <- gsva(
  gsvaParam(
  exprData = expr_matrix,      # normalized gene with gene symbols x sample matrix
  geneSets  = all_signatures,       # named list of gene vectors
  kcdf      = "Gaussian",      # "Gaussian" for log-CPM/TPM, "Poisson" for counts
  minSize   = 5,              # min genes per set
  maxSize   = 500              # max genes per set
  )
)


cat("GSVA scores dimensions:", dim(gsva_scores), "\n") #19 1013
head(gsva_scores)[,1:5]
rownames(gsva_scores) # gene set names
saveRDS(gsva_scores, "data/processed/immune/gsva_scores.rds")


## 4b. Single sample GSEA (ssGSEA) for immune deconvolution (Scores absolute, not relative to other samples);
##    How active is this pathway in this sample, regardless of other samples?

cat("\nRunning ssGSEA...\n")

ssgsea_scores <- gsva(
  ssgseaParam(
    exprData   = expr_matrix,
    geneSets   = all_signatures,
    minSize    = 5,
    maxSize    = 500
  )
)
cat("ssGSEA scores dimensions:", dim(ssgsea_scores), "\n")#19 1013 
saveRDS(ssgsea_scores, "data/processed/immune/ssgsea_scores.rds")



## ------------------------------------------------------------
## 5. Immune deconvolution using 7 different methods 
## 7 methods: quanTIseq,MCP-counter, EPIC,xCell, CIBERSORT, Estimate and Custom)
## ---------------------------------------------------------------------

# convert Ensembl IDs in the TPM matrix to gene symbols using the same mapping as above
rownames(tpm) <- sub("\\..*", "", rownames(tpm))
head(rownames(tpm))


tpm_filtered <- tpm[rownames(tpm) %in% id_map_1$ENSEMBL, ] # Subset TPM to only the Ensembl IDs that mapped successfully
cat("TPM dimensions after filtering:", dim(tpm_filtered), "\n") # 18568 1013

tpm_filtered <- tpm_filtered[id_map_1$ENSEMBL, ] #Reorder to match id_map_1 exactly
rownames(tpm_filtered) <- id_map_1$SYMBOL #Assign gene symbols
head(rownames(tpm_filtered))