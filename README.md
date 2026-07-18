# A Six-Gene Lipid Metabolism Score Stratifies Prognosis in Dedifferentiated Liposarcoma

## Overview
This repository contains all analysis code and processed data for the development and validation of a lipid metabolism-based prognostic score in dedifferentiated liposarcoma (DDLPS).

## Data Sources

| Dataset | n | Source | Accession |
|---------|---|--------|-----------|
| TCGA-SARC | 50 | The Cancer Genome Atlas | https://portal.gdc.cancer.gov/ |
| GSE21122 | 158 | Gene Expression Omnibus | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21122 |
| GSE159659 | 45 | Gene Expression Omnibus | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159659 |

## Core Genes
The six-gene lipid metabolism score comprises: **FASN**, **ACLY**, **SCD**, **SREBF1**, **HMGCR**, and **PPARG**.

## Repository Structure

```
├── scripts/              # All analysis scripts (Figures 1-6 + Supplementary)
├── results/
│   ├── figures/          # Main and supplementary figures (PDF/PNG/TIFF)
│   └── tables/           # Supplementary tables (S1-S4) + Key Resources Table
├── raw/                  # Processed expression matrices and metadata
├── data/                 # Final gene lists and pathway results
├── sessionInfo.txt       # R session information
├── .gitignore
└── README.md
```

## Requirements

- **R** >= 4.2.0
- R packages are documented in `sessionInfo.txt`

## Quick Start

```r
# Set working directory
setwd("path/to/Liposarcoma/2")

# Run scripts in order
source("scripts/Figure1_script.R")   # Study overview and cohort description
source("scripts/Figure2_script.R")   # Six-gene score development
source("scripts/Figure3_script.R")   # Survival and ROC analyses
source("scripts/Figure4_script.R")   # Pathway enrichment and immune analysis
source("scripts/Figure5_script.R")   # Nomogram and calibration
source("scripts/Figure6_script.R")   # External validation
```

## Supplementary Figures

```r
source("scripts/Supplementary_Figure_S1.R")
source("scripts/Supplementary_Figure_S2.R")
source("scripts/Supplementary_Figure_S3.R")
source("scripts/Supplementary_Figure_S4-S7.R")
```

## License

MIT License

## Citation

If you use this code or data, please cite: [Manuscript citation - to be added upon publication]
