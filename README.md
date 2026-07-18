# A Six-Gene Lipid Metabolism Score for DDLPS Prognosis

## Study Design

This study follows a discovery-validation-mechanism-translation framework:

- **Stage 1**: Data preparation and cohort overview
- **Stage 2**: Prognostic model development and validation
- **Stage 3**: Biological mechanism (PPARG/CD36 pathway, immune landscape)
- **Stage 4**: Clinical translation (nomogram, external validation)

## Key Questions Addressed

| Question | Script | Key Result |
|----------|--------|------------|
| Can the score independently predict prognosis? | stage2_prognostic_model | HR=5.83, P<0.01 (male) |
| What is the biological basis of the score? | stage3_mechanism | PPARG/CD36 pathway activation |
| How does the score affect immune microenvironment? | stage3_mechanism | Immune suppression |
| Is the score stable across independent cohorts? | stage4_clinical_translation | 3-cohort validation AUC>0.65 |

## Quick Start

```r
source("MAIN_PIPELINE.R")
```

## Repository Structure

- scripts/ — Analysis scripts organized by scientific stage
- results/figures/ — Main and supplementary figures
- results/tables/ — Supplementary tables
- raw/ — Processed expression matrices

## Data Sources

| Dataset | n | Source |
|---------|---|--------|
| TCGA-SARC | 50 | https://portal.gdc.cancer.gov/ |
| GSE21122 | 158 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21122 |
| GSE159659 | 45 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159659 |

## Citation

[To be updated upon publication]
