############################################################
# Normalization, filtering, and quality control 
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Define Immune signatures from MSigDB and literature
#   - Convert Ensembl IDs to gene symbols
#   - Run gene set variation analysis (GSVA) and single sample GSEA (ssGSEA) for immune signatures
#   - Run immune deconvolution using quanTIseq, MCP-counter, EPIC, xCell, and CIBERSORT 
#   - Combine all resutls for further analysis and visualization in downstream scripts
#   - Created heat map with key immune features



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
    library(tidyverse)
    library(RColorBrewer)
})


#create output directory
dir.create("results/figures/immune", recursive=TRUE, showWarnings = FALSE)
dir.create("data/processed/immune", recursive=TRUE, showWarnings = FALSE)
dir.create("data/immune", recursive=TRUE, showWarnings = FALSE)


## -----------------------------
## 1. Load Data
## -----------------------------
expr <- readRDS("data/processed/expression_vst_normalized.rds") #  DESeq2 VST as primary normalized count data 
cat("Expression data loaded. Dimensions:", dim(expr), "\n") #23059 1013

tpm <- readRDS("data/processed/tpm_filtered.rds")#23059 1013 
cat("TPM dimensions:", dim(tpm), "\n")

clinical <- readRDS("data/processed/clinical_filtered.rds")#1013 19
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
# Strip version suffix from rownames
rownames(expr) <- sub("\\..*", "", rownames(expr))
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

dup_genes<- rownames(expr_filtered_symbol)[duplicated(rownames(expr_filtered_symbol))]
cat("Number of duplicated gene symbols:", length(dup_genes), "\n") #0

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
S1<- data.frame(
  gene_set       = names(effective_sizes),
  original_size  = sapply(all_signatures, length),
  effective_size = effective_sizes
) 
print(S1)



## -----------------------------
## 4. GENE SET ENRICHMENT ANALYSIS: 
## Each sample gets a score for each gene set based on its expression of the genes in that set
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
saveRDS(gsva_scores, "results/immune/gsva_scores.rds")

#------------------------------------------------------------------------------
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
saveRDS(ssgsea_scores, "results/immune/ssgsea_scores.rds")



## ------------------------------------------------------------
## 5. Immune deconvolution using 5 different methods 
## 5 methods: quanTIseq,MCP-counter, EPIC,xCell, and CIBERSORT)
## ---------------------------------------------------------------------

# convert Ensembl IDs in the TPM  to gene symbols using the same mapping as above
rownames(tpm) <- sub("\\..*", "", rownames(tpm))
head(rownames(tpm))
dim(tpm)#23059  1013

tpm_filtered_symbol <- tpm[id_map_1$ENSEMBL, ] # Subset TPM to only the Ensembl IDs that mapped successfully; same rows, same order
rownames(tpm_filtered_symbol) <- id_map_1$SYMBOL        # same gene symbols as expr_filtered_symbol


cat("expr dimensions:", dim(expr_filtered_symbol), "\n")  # 18568 1013
cat("TPM dimensions: ", dim(tpm_filtered_symbol),  "\n")  # 18568 1013 

# Verify they are perfectly aligned
identical(rownames(expr_filtered_symbol), rownames(tpm_filtered_symbol))  # TRUE
identical(colnames(expr_filtered_symbol), colnames(tpm_filtered_symbol))  # TRUE


saveRDS(tpm_filtered_symbol, "data/processed/immune/tpm_filtered_symbol.rds")

#----------------------------------------------------------------------------------
# Method 1: quanTIseq
cat("\nRunning quanTIseq deconvolution...\n")

quantiseq_results <- deconvolute_quantiseq(
  gene_expression = tpm_filtered_symbol, #  genes x samples, linear TPM
  tumor = TRUE, # activates tumor-mode deconvolution
  arrays = FALSE, #RNA-seq data, not microarray
  scale_mrna = TRUE #corrects for cell-type mRNA content differences
)
cat("quantiseq_results dimensions: ", dim(quantiseq_results), "\n") #  11 1013
rownames(quantiseq_results)
head(quantiseq_results)[1:5,1:5]
saveRDS(quantiseq_results, "results/immune/quantiseq_results.rds")

# Summarize results across samples for each cell type
quantiseq_summary <- as.data.frame(t(quantiseq_results)) %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "fraction") %>%
  group_by(cell_type) %>%
  summarise(
    n        = n(),
    mean     = round(mean(fraction, na.rm = TRUE), 4),
    sd       = round(sd(fraction, na.rm = TRUE), 4),
    median   = round(median(fraction, na.rm = TRUE), 4),
    IQR      = round(IQR(fraction, na.rm = TRUE), 4),
    min      = round(min(fraction, na.rm = TRUE), 4),
    max      = round(max(fraction, na.rm = TRUE), 4),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean))
quantiseq_summary
write.csv(quantiseq_summary, "results/immune/quantiseq_summary.csv", row.names = FALSE)

#-----------------------------------------------------------------------------
# Method 2: MCP-counter:Scores are relative, not absolute
cat("\nRunning MCP-counter...\n")
mcp_results <- deconvolute_mcp_counter(
  gene_expression =tpm_filtered_symbol # genes x samples, linear TPM
)
cat("Dimensions of MCP-counter results: ", dim(mcp_results), "\n") # 10 1013
rownames(mcp_results)
head(mcp_results)[1:5,1:5]
saveRDS(mcp_results, "results/immune/mcp_results.rds")

#summarize results across samples for each cell type
mcp_summary <- as.data.frame(t(mcp_results)) %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "score") %>%
  group_by(cell_type) %>%
  summarise(
    n        = n(),
    mean     = round(mean(score, na.rm = TRUE), 4),
    sd       = round(sd(score, na.rm = TRUE), 4),
    median   = round(median(score, na.rm = TRUE), 4),
    IQR      = round(IQR(score, na.rm = TRUE), 4),
    min      = round(min(score, na.rm = TRUE), 4),
    max      = round(max(score, na.rm = TRUE), 4),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean))
mcp_summary
write.csv(mcp_summary, "results/immune/mcp_summary.csv", row.names = FALSE)

#-----------------------------------------------------------------------
# Method 3: EPIC:Fractions are absolute and sum to 1, but only for the cell types in the reference (B cells, CAFs, endothelial cells, macrophages, CD4 T cells, CD8 T cells, NK cells, and other cells)
cat("\nRunning EPIC...\n")
epic_results <- deconvolute_epic(
  gene_expression = tpm_filtered_symbol, ## genes x samples, linear TPM
  tumor = TRUE, ## use tumor-specific reference profiles
  scale_mrna = TRUE
)
cat("Dimensions of EPIC results: ", dim(epic_results), "\n") # 8 1013
rownames(epic_results)
head(epic_results)[1:5,1:5]
saveRDS(epic_results, "results/immune/epic_results.rds")

#summarize results across samples for each cell type
epic_summary <- as.data.frame(t(epic_results)) %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "fraction") %>%
  group_by(cell_type) %>%
  summarise(
    n        = n(),
    mean     = round(mean(fraction, na.rm = TRUE), 4),
    sd       = round(sd(fraction, na.rm = TRUE), 4),
    median   = round(median(fraction, na.rm = TRUE), 4),
    IQR      = round(IQR(fraction, na.rm = TRUE), 4),
    min      = round(min(fraction, na.rm = TRUE), 4),
    max      = round(max(fraction, na.rm = TRUE), 4),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean))
epic_summary
write.csv(epic_summary, "results/immune/epic_summary.csv", row.names = FALSE)


#---------------------------------------------------------------------------
# Method 4: xCell:Scores are enrichment-based, not fractions, and are relative to other samples, but cover a very wide range of cell types (67   immune and stromal cell types)
cat("\nRunning xCell...\n")
xcell_results <- deconvolute_xcell(
  gene_expression = tpm_filtered_symbol, ## genes x samples, linear TPM
  arrays = FALSE #RNA-seq mode
)
cat("Dimensions of xCell results: ", dim(xcell_results), "\n") # 67 1013
rownames(xcell_results)
head(xcell_results)[1:5,1:5]
saveRDS(xcell_results, "results/immune/xcell_results.rds")


xcell_results<- readRDS("results/immune/xcell_results.rds")
dim(xcell_results) # 67 1013
head(xcell_results)[,1:5]



#sumarize results across samples for each cell type
xcell_summary <- as.data.frame(t(xcell_results)) %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "score") %>%
  group_by(cell_type) %>% 
  summarise(
    n        = n(),
    mean     = round(mean(score, na.rm = TRUE), 4),
    sd       = round(sd(score, na.rm = TRUE), 4),
    median   = round(median(score, na.rm = TRUE), 4),
    IQR      = round(IQR(score, na.rm = TRUE), 4),
    min      = round(min(score, na.rm = TRUE), 4),
    max      = round(max(score, na.rm = TRUE), 4),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean))
xcell_summary
tail(xcell_summary, n=10)
write.csv(xcell_summary, "results/immune/xcell_summary.csv", row.names = FALSE)

#Summary of scores
summary_scores <- xcell_results[c("ImmuneScore",
                                   "StromaScore",
                                   "MicroenvironmentScore"), ]

xcell_score_summary <-as.data.frame(t(summary_scores)) %>%
  pivot_longer(everything(), names_to = "score_type", values_to = "score") %>%
  group_by(score_type) %>%
  summarise(
    mean   = round(mean(score,   na.rm = TRUE), 4),
    sd     = round(sd(score,     na.rm = TRUE), 4),
    median = round(median(score, na.rm = TRUE), 4),
    min    = round(min(score,    na.rm = TRUE), 4),
    max    = round(max(score,    na.rm = TRUE), 4),
    .groups = "drop"
  )
xcell_score_summary
write.csv(xcell_score_summary, "results/immune/xcell_score_summary.csv", row.names = FALSE)

#---------------------------------------------------------------------------
#Method 5: CIBERSORT: Fractions are relative to other samples and only cover 22 immune cell types
##Permutation=1000
cat("\nRunning CIBERSORT...\n")

#Filter to CIBERSORT LM22 signature genes only
# Download LM22.txt from CIBERSORTx portal
                   
lm22 <- read.table("data/immune/LM22.txt",sep= "\t",header = TRUE, row.names = 1,check.names = FALSE)
cat("Dimension of lm22:", dim(lm22), "\n") # 547 22             
                   
#Extract overlapping genes                  
lm22_genes    <- rownames(lm22)
overlap_genes <- intersect(lm22_genes, rownames(tpm_filtered_symbol))

cat("LM22 genes found in your matrix:", length(overlap_genes), "\n")  # 435 

tpm_cibersort_overlap <- tpm_filtered_symbol[overlap_genes, ]
cat("Final dimensions:", dim(tpm_cibersort_overlap), "\n")  # 435  x 1013
head(tpm_cibersort_overlap)[1:5,1:5]

# tpm_cibersort <-as.data.frame(tpm_cibersort)
# tpm_cibersort <- tibble::rownames_to_column(tpm_cibersort, "GeneSymbol")
# write.table(tpm_cibersort,"data/processed/immune/tpm_cibersort.txt",sep= "\t",row.names = FALSE, quote = FALSE)

#checking missing genes 
missing <- setdiff(lm22_genes, rownames(tpm_filtered_symbol))
cat("Missing genes:", length(missing), "\n")
head(missing,20)

#case mismatch check
missing_upper <- toupper(missing)
matrix_upper  <- toupper(rownames(tpm_filtered_symbol))
case_matches  <- sum(missing_upper %in% matrix_upper)
cat("Recoverable by case fix:", case_matches, "\n")

lost_in_mapping <- setdiff(missing,
                            rownames(tpm_filtered_symbol))
cat("Genes lost in ID mapping:", length(lost_in_mapping), "\n")

dup_genes_cibersort <- rownames(tpm_filtered_symbol)[duplicated(rownames(tpm_filtered_symbol))]
cat("Number of duplicated gene symbols in CIBERSORT matrix:", length(dup_genes_cibersort), "\n") #0


# 1013 samples split into 5 batches
batch_plan <- list(
  batch_1 = 1:200,
  batch_2 = 201:400,
  batch_3 = 401:600,
  batch_4 = 601:800,
  batch_5 = 801:1013   # 213 samples
)


dir.create("data/processed/cibersort", 
           recursive = TRUE, showWarnings = FALSE)

for (batch_name in names(batch_plan)) {
  idx      <- batch_plan[[batch_name]]
  batch_df <- tpm_cibersort_overlap[, idx] %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneSymbol")


write.table(
    batch_df,
    sprintf("data/processed/cibersort/%s.txt", batch_name),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  cat(sprintf("%s: samples %4d - %4d (%d samples)\n",
              batch_name, min(idx), max(idx), length(idx)))
}

# Combina all the downloaded CIBERSORTx result CSVs

result_files <- list.files(
  "data/processed/cibersort/results",
  pattern    = "\\.csv$",
  full.names = TRUE
)

cat("Result files found:", length(result_files), "\n")  # should be 6
head(result_files)


cibersort_merged <- result_files %>% 
  lapply(read.csv, check.names = FALSE) %>% 
  bind_rows()

cat("Total samples merged:", nrow(cibersort_merged), "\n")  # should be 1013

# Verify no samples duplicated
cat("Unique samples:", n_distinct(cibersort_merged[, 1]), "\n")  # should also be 1013
cat("Dimensions of merged CIBERSORT results:", dim(cibersort_merged), "\n")  # should be 1013 x 26 (1 barcode + 22 cell types)
saveRDS(cibersort_merged,"results/immune/cibersort_results_merged.rds")
head(cibersort_merged)[,1:5]
colnames(cibersort_merged)

# Verify sample order matches your metadata
dim(cibersort_merged)#1013   26
cibersort_merged_1 <- cibersort_merged %>% 
  rename(barcode = "Mixture") %>% 
  arrange(match(barcode, colnames(tpm_filtered_symbol)))  # reorder to match TPM
head(cibersort_merged_1)[,1:5]
head(tpm_filtered_symbol)[,1:5]

# Save merged results
saveRDS(cibersort_merged_1,"results/immune/cibersort_results_merged_barcode.rds")
write.csv(cibersort_merged_1,"data/processed/immune/cibersort_results_merged_barcode.csv", row.names = FALSE)




## ------------------------------------------------------------
## 6. Combine and compare results across methods
## ---------------------------------------------------------------------
# Transpose deconvolution results to have samples as rows

#gsva
gsva_scores  <- readRDS("results/immune/gsva_scores.rds")
dim(gsva_scores) # 19 1013
head(gsva_scores)[,1:5]
rownames(gsva_scores) # gene set names in rows, samples in columns

gsva_df <- as.data.frame(t(gsva_scores))
head(gsva_df)[,1:5]

colnames(gsva_df) <- paste0("GSVA_", colnames(gsva_df))
head(gsva_df)[,1:5]

#ssgsea
ssgsea_scores <- readRDS("results/immune/ssgsea_scores.rds")
dim(ssgsea_scores) # 19 1013
head(ssgsea_scores)[,1:5]
rownames(ssgsea_scores) # gene set names in rows, samples in columns

ssgsea_df <- as.data.frame(t(ssgsea_scores))
head(ssgsea_df)[,1:5]
colnames(ssgsea_df) <- paste0("ssGSEA_", colnames(ssgsea_df))
head(ssgsea_df)[,1:5]

#quanTIseq
quantiseq_results <- readRDS("results/immune/quantiseq_results.rds")
dim(quantiseq_results) # 11 1013
head(quantiseq_results)[,1:5]
rownames(quantiseq_results) # cell types in rows, samples in columns

quantiseq_df <- as.data.frame(t(quantiseq_results))
head(quantiseq_df)[,1:5]
colnames(quantiseq_df) <- paste0("quanTIseq_", colnames(quantiseq_df))
head(quantiseq_df)[,1:5]

#MCP-counter
mcp_results <- readRDS("results/immune/mcp_results.rds")
dim(mcp_results) # 10 1013
head(mcp_results)[,1:5]
rownames(mcp_results) # cell types in rows, samples in columns

mcp_df <- as.data.frame(t(mcp_results))
colnames(mcp_df) <- paste0("MCP_", colnames(mcp_df))
head(mcp_df)[,1:5]

#EPIC
epic_results <- readRDS("results/immune/epic_results.rds")
dim(epic_results) # 8 1013
head(epic_results)[,1:5]
rownames(epic_results) # cell types in rows, samples in columns

epic_df <- as.data.frame(t(epic_results))
colnames(epic_df) <- paste0("EPIC_", colnames(epic_df))
head(epic_df)[,1:5]
dim(epic_df) # 1013 8

#Xcell
xcell_results <- readRDS("results/immune/xcell_results.rds")
dim(xcell_results) # 67 1013
head(xcell_results)[,1:5]
rownames(xcell_results) # cell types in rows, samples in columns

xcell_df <- as.data.frame(t(xcell_results))
colnames(xcell_df) <- paste0("xCell_", colnames(xcell_df))
head(xcell_df)[,1:5]
dim(xcell_df) # 1013 67



#cibersort
cibersort_results <- readRDS("results/immune/cibersort_results_merged_barcode.rds")
dim(cibersort_results) # 1013 26 
head(cibersort_results)[,1:5]


cibersort_df <- cibersort_results %>%
  column_to_rownames("barcode") %>%
  as.data.frame()
head(cibersort_df)[,1:5]

colnames(cibersort_df) <- paste0("CIBERSORT_", colnames(cibersort_df))
head(cibersort_df)[,1:5]
dim(cibersort_df) # 1013  25 (1 barcode + 22 cell types + 2 summary scores)



# Combine all features
immune_features <- cbind(
  gsva_df,
  ssgsea_df,
  quantiseq_df,
  mcp_df,
  epic_df,
  xcell_df,
  cibersort_df
)
dim(immune_features) # 1013  159  (19+19+11+10+8+67+22) = 156 features
head(immune_features)[,1:5]


# Ensure sample IDs match
immune_features <- cbind(sample = colnames(expr), immune_features)

saveRDS(immune_features, "data/processed/immune/immune_features_combined.rds")

cat("\nTotal immune features:", ncol(immune_features), "\n")

colnames(immune_features)[c(1:10,159:160)]




# Correlation heatmap of key immune features from different methods
key_features <- c(
  # CD8 T cells (cytotoxic)
  "GSVA_T_cells_CD8", "ssGSEA_T_cells_CD8", "quanTIseq_T.cells.CD8",
  "MCP_CD8 T cells", "EPIC_CD8_Tcells", "CIBERSORT_T cells CD8",

  # Tregs (immunosuppressive)
  "GSVA_T_reg", "ssGSEA_T_reg", "quanTIseq_Tregs",
  "CIBERSORT_T cells regulatory (Tregs)", "xCell_Tregs",

  # T cell exhaustion
  "GSVA_T_exhaustion", "ssGSEA_T_exhaustion",

  # B cells
  "GSVA_B_cells", "ssGSEA_B_cells", "quanTIseq_B.cells",
  "MCP_B lineage", "EPIC_Bcells", "CIBERSORT_B cells memory",

  # Macrophages M1 vs M2
  "GSVA_Macrophages_M1", "GSVA_Macrophages_M2",
  "quanTIseq_Macrophages.M1", "quanTIseq_Macrophages.M2",
  "EPIC_Macrophages", "CIBERSORT_Macrophages M1", "CIBERSORT_Macrophages M2",
  "xCell_Macrophages M1", "xCell_Macrophages M2",

  # NK cells
  "quanTIseq_NK.cells", "MCP_NK cells", "EPIC_NKcells",

  # Dendritic cells
  "GSVA_Dendritic_cells", "ssGSEA_Dendritic_cells",
  "quanTIseq_Dendritic.cells", "MCP_Myeloid dendritic cells",

  # Inflammatory signatures
  "GSVA_HALLMARK_INFLAMMATORY_RESPONSE", "ssGSEA_HALLMARK_INFLAMMATORY_RESPONSE",
  "GSVA_HALLMARK_INTERFERON_GAMMA_RESPONSE", "ssGSEA_HALLMARK_INTERFERON_GAMMA_RESPONSE",

  # Cytotoxicity & immunosuppression
  "GSVA_Cytotoxicity", "ssGSEA_Cytotoxicity",
  "GSVA_Immunosuppression", "ssGSEA_Immunosuppression",

  # Stromal
  "MCP_Endothelial cells", "MCP_Fibroblasts",
  "EPIC_CAFs", "EPIC_Endothelial",
  "xCell_ImmuneScore", "xCell_StromaScore", "xCell_MicroenvironmentScore"
)


cor_matrix <- immune_features %>%
  dplyr::select(all_of(key_features)) %>%
  as.matrix()

dim(cor_matrix)#1013   50

#Calculate spearman coorelation:CIBERSORT and quanTIseq outputs are proportions (0–1) while GSVA/ssGSEA scores are continuous enrichment scores.
#Spearman handles this mixed scale better and is more robust to the zero-inflation common in deconvolution outputs.
cor_result <- cor(cor_matrix, method = "spearman", use = "pairwise.complete.obs")
dim(cor_result)#50 50             

png("results/figures/immune/immune_features_correlation.png", 
    width = 14, height = 14, units = "in", res = 300)
pheatmap(
  cor_result,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-1, 1, length.out = 101),
  clustering_method = "ward.D2",
  treeheight_row    = 40,
  treeheight_col    = 40,
  fontsize          = 7,
  fontsize_row      = 6,
  fontsize_col      = 6,
  main              = "Spearman Correlation of Immune Features",
  border_color      = NA,
  display_numbers   = FALSE
)
dev.off()









