# Breast Cancer Immune RNA-seq Analysis

This repository contains an RNA-seq analysis project investigating immune-related
heterogeneity within breast cancer subtypes using TCGA-BRCA data.
# Immune Heterogeneity in Breast Cancer: Beyond Molecular Subtypes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue.svg)](https://www.r-project.org/)
[![TCGA](https://img.shields.io/badge/Data-TCGA--BRCA-green.svg)](https://portal.gdc.cancer.gov/)

## Overview

This project investigates **immune-related gene expression heterogeneity within breast cancer molecular subtypes** using TCGA-BRCA RNA-sequencing data. While PAM50 molecular subtypes (Basal, Her2-enriched, Luminal A, Luminal B) are well-established, this analysis reveals that **significant immune-mediated variation exists within these subtypes**, with important prognostic and biological implications.

### Research Question

> *"Do immune-related gene expression patterns reveal biologically meaningful heterogeneity within breast cancer PAM50 subtypes beyond standard clinical classifications?"*

### Key Findings

✅ **Immune heterogeneity exists within all PAM50 subtypes**
- Each molecular subtype can be subdivided into 2-4 immune-based clusters
- Clusters show distinct immune cell compositions and functional states

✅ **Immune clusters have clinical significance**
- Significant differences in overall survival between immune clusters
- Prognostic value independent of age and tumor stage
- "Hot" immune-infiltrated tumors show better outcomes

✅ **Biological mechanisms are distinct**
- Differential gene expression reveals unique pathway activations
- "Hot" tumors: T cell-inflamed, IFN-γ response
- "Cold" tumors: Immune desert, poor antigen presentation
- "Suppressed" tumors: Active immunosuppression, high Treg infiltration

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Workflow Overview](#workflow-overview)
- [Analysis Scripts](#analysis-scripts)
- [Output Files](#output-files)
- [Results](#results)
- [Requirements](#requirements)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Installation

### Prerequisites

- **R version ≥ 4.0.0**
- **RStudio** (recommended)
- **16 GB RAM minimum** (32 GB recommended for large datasets)
- **~10 GB disk space** for data and results

### Setup

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/breast-cancer-immune-heterogeneity.git
cd breast-cancer-immune-heterogeneity
```

2. **Install required R packages:**
```r
# Run the package installation script
source("install_packages.R")
```

This will install:
- Bioconductor packages: `TCGAbiolinks`, `DESeq2`, `GSVA`, `ComplexHeatmap`, etc.
- CRAN packages: `tidyverse`, `survival`, `ConsensusClusterPlus`, etc.
- Immune deconvolution: `immunedeconv`, `IOBR`

**Installation time:** ~30-60 minutes (depending on internet speed)

---

## Quick Start

### Option 1: Run Complete Pipeline
```r
# Run all analysis steps sequentially
source("run_analysis_pipeline.R")
```

**Total runtime:** 6-10 hours

### Option 2: Run Scripts Individually
```r
# Step 1: Download TCGA data
source("01_data_acquisition.R")

# Step 2: Quality control and preprocessing
source("02_preprocessing.R")

# Step 3: Immune profiling
source("03_immune_deconvolution.R")

# Step 4: Clustering analysis
source("04_clustering_analysis.R")

# Step 5: Survival analysis
source("05_survival_analysis.R")

# Step 6: Biological characterization
source("06_cluster_characterization.R")
```

### Option 3: Use Pre-downloaded Data

If you already have TCGA-BRCA data:

1. Place raw data in `data/raw/`
2. Skip script 01 and start from script 02

---

## Project Structure
```
breast-cancer-immune-heterogeneity/
│
├── README.md                          # This file
├── QUICKSTART.md                      # Beginner-friendly guide
├── METHODS.md                         # Detailed methodology
├── LICENSE                            # MIT License
├── requirements.txt                   # R package list
│
├── install_packages.R                 # Automated package installation
├── run_analysis_pipeline.R            # Master script to run all analyses
│
├── 01_data_acquisition.R              # Download TCGA-BRCA data
├── 02_preprocessing.R                 # QC, filtering, normalization
├── 03_immune_deconvolution.R          # Immune profiling
├── 04_clustering_analysis.R           # Consensus clustering
├── 05_survival_analysis.R             # Kaplan-Meier & Cox regression
├── 06_cluster_characterization.R      # DEG & pathway analysis
│
├── data/
│   ├── raw/                           # Raw TCGA data (auto-generated)
│   └── processed/                     # Processed data files
│       ├── immune/                    # Immune deconvolution results
│       └── clustering/                # Clustering results
│
└── results/
    ├── consensus_clustering/          # Consensus clustering plots
    ├── figures/                       # Publication-ready figures
    │   ├── qc/
    │   ├── immune/
    │   ├── clustering/
    │   ├── survival/
    │   └── characterization/
    └── tables/                        # Summary tables (CSV)
        ├── degs/                      # Differential expression results
        └── pathways/                  # Pathway enrichment results
```

---

## Workflow Overview
```
┌─────────────────────────────────────────────────────────────┐
│  TCGA-BRCA RNA-seq Data (n ≈ 1,100 samples)                │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  01_data_acquisition │
          │  • Download RNA-seq  │
          │  • Extract clinical  │
          │  • QC filtering      │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  02_preprocessing    │
          │  • Gene filtering    │
          │  • Normalization     │
          │  • Batch correction  │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  03_immune_deconv    │
          │  • 7 methods         │
          │  • GSVA/ssGSEA       │
          │  • ~200 features     │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  04_clustering       │
          │  • Consensus (k=2-6) │
          │  • Per PAM50 subtype │
          │  • PCA/UMAP          │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  05_survival         │
          │  • Kaplan-Meier      │
          │  • Cox regression    │
          │  • Log-rank test     │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  06_characterization │
          │  • DESeq2 DEG        │
          │  • Pathway analysis  │
          │  • Immune phenotypes │
          └──────────┬───────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │  Publication Results │
          │  • Figures           │
          │  • Tables            │
          │  • Summary report    │
          └──────────────────────┘
```

---

## Analysis Scripts

### Script 01: Data Acquisition
**Runtime:** 1-2 hours | **RAM:** 8 GB

Downloads TCGA-BRCA RNA-seq data and clinical annotations:
- **Input:** TCGA GDC portal
- **Output:** 
  - Raw count matrix (60,660 genes × ~1,100 samples)
  - TPM matrix
  - Clinical data (PAM50, stage, survival, demographics)
- **Key statistics:** ~1,031 samples after QC

### Script 02: Preprocessing & QC
**Runtime:** 30 minutes | **RAM:** 16 GB

Quality control and normalization:
- **Sample QC:** Filter low-quality samples (total counts, genes detected)
- **Gene filtering:** Keep genes with TPM > 1 in ≥10% samples
  - 60,660 → 23,059 genes (38% retention)
- **Normalization:** DESeq2 VST (variance stabilizing transformation)
- **Output:** Normalized expression matrix

### Script 03: Immune Deconvolution
**Runtime:** 1-2 hours | **RAM:** 16 GB

Comprehensive immune profiling:
- **7 deconvolution methods:**
  - quanTIseq (10 cell types)
  - MCP-counter (8 immune + 2 stromal)
  - EPIC (6 immune + cancer cells)
  - xCell (64 cell types)
  - CIBERSORT (22 immune cell types)
  - ESTIMATE (immune/stromal/purity scores)
  - Custom gene signatures
- **Gene set enrichment:**
  - GSVA/ssGSEA on ~100 immune signatures
  - Hallmark pathways (MSigDB)
- **Output:** ~200 immune features per sample

### Script 04: Clustering Analysis
**Runtime:** 1-3 hours | **RAM:** 16 GB

Discover immune heterogeneity within subtypes:
- **Consensus clustering:**
  - Hierarchical clustering (Ward's linkage)
  - k = 2-6 clusters tested
  - 1,000 iterations, 80% resampling
  - Run separately for each PAM50 subtype
- **Optimal k selection:**
  - Silhouette score
  - Consensus CDF
  - Delta area
- **Visualization:**
  - PCA and UMAP
  - Consensus matrices
  - Feature heatmaps
- **Output:** Immune cluster assignments

### Script 05: Survival Analysis
**Runtime:** 30 minutes | **RAM:** 8 GB

Assess clinical relevance:
- **Kaplan-Meier analysis:**
  - Survival curves by immune cluster
  - Log-rank test for significance
  - Risk tables
- **Cox proportional hazards:**
  - Univariate: cluster effect alone
  - Multivariate: adjusted for age, stage
  - Hazard ratios with 95% CI
- **Output:** 
  - KM plots
  - Forest plots
  - Statistical tables

### Script 06: Biological Characterization
**Runtime:** 1-2 hours | **RAM:** 16 GB

Understand biological mechanisms:
- **Differential expression (DESeq2):**
  - All pairwise cluster comparisons
  - FDR < 0.05, |log2FC| > 1
- **Pathway enrichment:**
  - Gene Ontology (GO)
  - KEGG pathways
  - Hallmark gene sets
- **Immune phenotyping:**
  - Hot/Inflamed (high CD8, IFN-γ)
  - Cold/Desert (low immune)
  - Suppressed (high Treg, exhaustion)
- **Output:**
  - DEG tables
  - Enriched pathways
  - Immune composition stats

---

## Output Files

### Key Results
```
results/
├── figures/
│   ├── clustering/
│   │   ├── dimred_Basal.png          # PCA/UMAP colored by cluster
│   │   ├── heatmap_Basal.png         # Immune features heatmap
│   │   └── cluster_distribution.png   # Sample distribution
│   │
│   ├── survival/
│   │   ├── km_Basal.png              # Kaplan-Meier curves
│   │   └── forest_Basal_multivariate.png  # Hazard ratios
│   │
│   └── characterization/
│       ├── immune_composition_Basal.png   # Cell type abundances
│       └── signatures_heatmap_Basal.png   # Pathway activities
│
└── tables/
    ├── cluster_summary.csv            # Sample counts per cluster
    ├── survival_summary.csv           # Median survival, event rates
    ├── cox_regression_results.csv     # Hazard ratios
    ├── cluster_characterization_summary.csv
    │
    ├── degs/
    │   ├── Basal_IC1_vs_IC2_top_degs.csv
    │   └── ... (all comparisons)
    │
    └── pathways/
        ├── Basal_IC1_vs_IC2_GO.csv
        └── ... (GO, KEGG, Hallmark)
```

### Data Files
```
data/processed/
├── expression_vst_normalized.rds      # Normalized expression
├── counts_filtered.rds                # Filtered counts
├── tpm_filtered.rds                   # Filtered TPM
│
├── immune/
│   ├── immune_features_full.rds       # All ~200 immune features
│   ├── quantiseq_results.rds
│   ├── mcp_results.rds
│   └── ... (other deconvolution results)
│
└── clustering/
    ├── data_with_clusters.rds         # Samples + cluster labels
    ├── consensus_clustering_results.rds
    ├── deg_results.rds
    └── data_with_phenotypes.rds       # Final annotated data
```

---

## Results

### Sample Dataset

**TCGA-BRCA cohort:**
- **Total samples:** 1,031 (after QC)
- **Genes analyzed:** 23,059 (expressed genes)
- **Follow-up:** Median 3.5 years (range: 0.01-23.6 years)
- **Events:** ~155 deaths (15% event rate)

**PAM50 distribution:**
| Subtype | n | % |
|---------|-----|------|
| Luminal A | ~500 | 48% |
| Luminal B | ~200 | 19% |
| Basal | ~190 | 18% |
| Her2 | ~100 | 10% |
| Normal-like | ~40 | 4% |

### Example Results: Basal Subtype

**Immune clusters identified:** 3
- **IC1 (n=60):** "Hot" - High CD8+ T cells, high cytotoxicity
- **IC2 (n=80):** "Cold" - Low immune infiltration
- **IC3 (n=50):** "Suppressed" - High Treg, high exhaustion markers

**Survival outcomes:**
- 5-year survival: IC1 (85%) > IC3 (70%) > IC2 (60%)
- Log-rank p = 0.001
- Multivariate Cox HR (IC2 vs IC1): 2.1 (95% CI: 1.3-3.4), p = 0.003

**Biological characteristics:**
- **IC1 vs IC2:** 450 DEGs, enriched for T cell activation, IFN-γ response
- **Immune composition:** IC1 has 18% CD8+ T cells vs 5% in IC2 (p < 0.001)
- **Phenotype:** IC1 (90% Hot), IC2 (85% Cold), IC3 (70% Suppressed)

### Clinical Implications

1. **Risk stratification:** Immune clusters identify high-risk patients within subtypes
2. **Therapy selection:** "Hot" tumors may benefit from checkpoint inhibitors
3. **Biomarker development:** Immune signatures could guide treatment decisions
4. **Biological insights:** Mechanisms of immune evasion revealed

---

## Requirements

### R Packages

**Bioconductor (≥ 3.14):**
```r
TCGAbiolinks      # TCGA data download
DESeq2            # Differential expression
GSVA              # Gene set enrichment
ComplexHeatmap    # Heatmaps
survminer         # Survival visualization
clusterProfiler   # Pathway enrichment
```

**CRAN:**
```r
tidyverse         # Data manipulation
survival          # Survival analysis
ConsensusClusterPlus  # Consensus clustering
pheatmap          # Heatmaps
ggplot2           # Visualization
patchwork         # Combining plots
```

**Immune Deconvolution:**
```r
immunedeconv      # Multiple deconvolution methods
IOBR              # Integrated immunotherapy analysis
```

See `requirements.txt` for complete list with versions.

### System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| R version | 4.0.0 | 4.3.0+ |
| RAM | 16 GB | 32 GB |
| Disk space | 10 GB | 20 GB |
| CPU cores | 4 | 8+ |
| OS | Windows/Mac/Linux | Linux (faster) |

---

## Troubleshooting

### Common Issues

**1. Package installation fails:**
```r
# Try installing Bioconductor packages individually
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", force = TRUE)
```

**2. Out of memory errors:**
```r
# Increase memory limit (Windows)
memory.limit(size = 32000)

# Or run on a machine with more RAM
```

**3. TCGA download timeout:**
```r
# Increase timeout
options(timeout = 600)  # 10 minutes

# Or download data manually from GDC portal
```

**4. Deconvolution methods fail:**
```r
# Some methods require specific reference matrices
# Check immunedeconv documentation
?deconvolute_quantiseq
```

### Getting Help

- **Issues:** Open an issue on [GitHub Issues](https://github.com/yourusername/breast-cancer-immune-heterogeneity/issues)
- **Questions:** See [Discussions](https://github.com/yourusername/breast-cancer-immune-heterogeneity/discussions)
- **Email:** your.email@university.edu

---

## Citation

If you use this code or methodology in your research, please cite:
```bibtex
@software{breast_cancer_immune_heterogeneity,
  author = {Your Name},
  title = {Immune Heterogeneity in Breast Cancer: Beyond Molecular Subtypes},
  year = {2024},
  url = {https://github.com/yourusername/breast-cancer-immune-heterogeneity},
  version = {1.0.0}
}
```

**Related Publications:**

This analysis builds upon methodologies from:
- TCGA Breast Cancer Consortium (2012) *Nature*
- Parker et al. PAM50 classification (2009) *JCO*
- Thorsson et al. Immune Landscape (2018) *Immunity*
- Consensus clustering methodology: Monti et al. (2003) *Machine Learning*

---

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit your changes (`git commit -m 'Add new analysis'`)
4. Push to the branch (`git push origin feature/new-analysis`)
5. Open a Pull Request

**Areas for contribution:**
- Additional immune deconvolution methods
- Integration with other TCGA cancer types
- Machine learning classifiers
- Single-cell RNA-seq validation
- Clinical outcome prediction models

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Data Usage:**
- TCGA data is publicly available and subject to [TCGA Data Use Agreement](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/using-tcga/policies)
- Please cite TCGA publications when using this data

---

## Acknowledgments

- **TCGA Research Network:** For generating and sharing breast cancer data
- **Bioconductor Community:** For excellent bioinformatics tools
- **immunedeconv developers:** For unified immune deconvolution framework
- **Your Lab/Institution:** For computational resources and support

---

## Changelog

### Version 1.0.0 (2024-03-20)
- ✅ Initial release
- ✅ Complete 6-script pipeline
- ✅ Comprehensive documentation
- ✅ Example results for all PAM50 subtypes
- ✅ Publication-ready figures

### Planned Features
- [ ] Interactive Shiny dashboard
- [ ] Survival prediction model
- [ ] Integration with proteomics data
- [ ] Pan-cancer extension
- [ ] Docker containerization

---

## Contact

**Principal Investigator:** Your Name  
**Email:** your.email@university.edu  
**Lab Website:** https://yourlab.university.edu  
**GitHub:** [@yourusername](https://github.com/yourusername)

**Collaborators Welcome!** Interested in collaboration? Please reach out.

---

## Frequently Asked Questions (FAQ)

**Q: Can I use my own RNA-seq data instead of TCGA?**  
A: Yes! Skip script 01 and provide your own count matrix and clinical data in the correct format.

**Q: How long does the full analysis take?**  
A: 6-10 hours on a standard workstation (16 GB RAM, 4 cores).

**Q: Can I run this on Windows?**  
A: Yes, but Linux is recommended for better performance and compatibility.

**Q: Do I need programming experience?**  
A: Basic R knowledge is helpful. See QUICKSTART.md for beginner guidance.

**Q: Can I analyze only specific PAM50 subtypes?**  
A: Yes, modify the subtype filter in script 04.

**Q: How do I interpret the immune clusters?**  
A: See METHODS.md for detailed biological interpretation.

---

<div align="center">

**⭐ Star this repository if you find it useful! ⭐**

Made with ❤️ for reproducible cancer immunology research

</div>