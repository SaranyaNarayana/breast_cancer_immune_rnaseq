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




## -----------------------------
## 3. Gene set variation analysis
## -----------------------------

cat("\nRunning GSVA...\n")

# Convert expression to matrix 
expr_matrix <- as.matrix(expr)
head(expr_matrix)[,1:5] # has Ensembl IDs

head(all_signatures[[1]]) #has gene symbols, need to convert to Ensembl IDs

#Convert Ensembl IDs in the  expr_matrix to gene symbols using AnnotationDbi and org.Hs.eg.db
install.packages(c("AnnotationDbi", "org.Hs.eg.db"))
library(AnnotationDbi)
library(org.Hs.eg.db)  # install via BiocManager::install("org.Hs.eg.db")

# Strip version suffix from rownames
rownames(expr_matrix) <- sub("\\..*", "", rownames(expr_matrix))


id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rownames(expr_matrix),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)
head(id_map)
cat("ID mapping dimensions:", dim(id_map), "\n") #23267 2 

# Remove unmapped or duplicate symbols
id_map <- id_map[!is.na(id_map$SYMBOL) & !duplicated(id_map$SYMBOL), ]

# Update matrix rownames
expr_matrix <- expr_matrix[id_map$ENSEMBL, ]
rownames(expr_matrix) <- id_map$SYMBOL













#library(biomaRt)
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
dim(mart)


# Convert gene symbols in your gene sets to Ensembl IDs
converted <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = unique(unlist(gene_sets)),
  mart       = mart
)

# Rebuild gene sets with Ensembl IDs
gene_sets_ensembl <- lapply(gene_sets, function(gs) {
  matched <- converted$ensembl_gene_id[converted$hgnc_symbol %in% gs]
  matched[matched != ""]
})












# Build the parameter object
params <- gsvaParam(
  exprData = expr_matrix,      # normalized gene x sample matrix
  geneSets  = all_signatures,       # named list of gene vectors
  kcdf      = "Gaussian",      # "Gaussian" for log-CPM/TPM, "Poisson" for counts
  minSize   = 10,              # min genes per set
  maxSize   = 500              # max genes per set
)

cat("Number ofremaing gene sets:", length(gsva(params) |> rownames()))

# Run GSVA
gsva_scores <- gsva(params)

saveRDS(gsva_scores, "data/processed/immune/gsva_scores.rds")