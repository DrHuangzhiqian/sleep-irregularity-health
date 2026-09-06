# Accelerometer-Derived Sleep Irregularity and Human Health

## Overview

This repository provides the analytical code and reproducible scripts underlying a comprehensive investigation of accelerometer-derived sleep irregularity and its relevance to human health and disease.

Using large-scale data from the UK Biobank, this study evaluates accelerometer-derived sleep irregularity in relation to multi-system disease outcomes, metabolomic and proteomic signatures, inflammatory markers, genetic architecture, potential causal pathways, and longitudinal sleep regularity trajectories.

The analytical workflow includes sleep phenotype derivation, prospective epidemiological modelling, omics association analyses, mediation analyses, genome-wide association analyses, Mendelian randomization, and longitudinal change pattern analyses.

## Code Description

The main scripts included in this repository are:

| File | Description |
|---|---|
| `GGIR.R` | Processes raw accelerometer data using GGIR. |
| `Cox_main_analysis.R` | Performs Cox regression analyses of sleep regularity with a range of disease outcomes. |
| `Cox_combinedgroup_duration_normal.R` | Performs Cox analyses of joint sleep phenotype groups across outcomes. |
| `linear_omics.R` | Performs linear association analyses of omics features with sleep irregularity. |
| `linear_combinedgroup_omics.R` | Performs linear analyses of joint sleep phenotype groups with omics traits. |
| `Lasso_met_scores_min.R` | Builds LASSO-derived metabolomic signature scores for sleep irregularity. |
| `Lasso_pro_scores_min.R` | Builds LASSO-derived proteomic signature scores for sleep irregularity. |
| `Cox_signature_Met_Outcomes.R` | Performs Cox analyses of metabolomic signature scores across outcomes. |
| `Cox_signature_Pro_Outcomes.R` | Performs Cox analyses of proteomic signature scores across outcomes. |
| `signature_metabolic_med.R` | Performs mediation analyses using metabolomic signature scores. |
| `signature_proteomic_med.R` | Performs mediation analyses using proteomic signature scores. |
| `single_Metab_parado.R` | Performs single-metabolite mediation analyses. |
| `single_Prote_parado.R` | Performs single-protein mediation analyses. |
| `single_INFLA_parado.R` | Performs single inflammatory marker mediation analyses. |
| `score_INFLA_med.R` | Performs INFLA-score mediation analyses across outcomes. |
| `MR analysis.R` | Performs Mendelian randomization analyses of sleep irregularity and health outcomes. |
| `longitudinal_change_patterns.R` | Performs analyses of longitudinal change patterns in sleep regularity. |
| `regenie_step1_template.sh` | Provides a template for REGENIE step 1 genome-wide association analysis. |
| `regenie_step2_template.sh` | Provides a template for REGENIE step 2 chromosome-wise association testing. |
| `PRSCS_template.sh` | Provides a template for PRS-CS-auto posterior effect estimation. |
| `prs_plink2_template.sh` | Provides a template for PLINK2-based polygenic risk score calculation. |

## System Requirements

Analyses were performed on the UK Biobank Research Analysis Platform using:

- R software, version 4.4.1
- PLINK2 for genome-wide association analyses
- Additional R packages for epidemiological, omics, mediation, and Mendelian randomization analyses

The main R packages used include packages for data processing, survival analysis, regression modelling, visualization, penalized regression, mediation analysis, and Mendelian randomization.

No non-standard hardware is required for running small-scale scripts or analyses using the demonstration dataset. Full-scale analyses using individual-level UK Biobank data, omics data, or genome-wide genetic data require access to the UK Biobank Research Analysis Platform or other high-performance computing resources.

## Installation

Clone the repository:

```bash
git clone https://github.com/DrHuangzhiqian/sleep-irregularity-health.git
cd sleep-irregularity-health
````

Install required R packages according to the package calls listed at the beginning of each script.

For common R dependencies, users may install:

```r
install.packages(c(
  "data.table",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "ggplot2",
  "survival",
  "broom",
  "purrr",
  "glmnet"
))
```

Specialized packages for Mendelian randomization or genetic analyses should be installed according to their official documentation.

Typical installation time for common R packages is approximately 10–30 minutes on a standard desktop computer.

## Reproducibility

This repository contains analytical code only and does not include any individual-level UK Biobank data, participant identifiers, omics data, genetic data, or other restricted participant-level records.

Researchers who wish to test the code using simulated data for workflow validation may contact us by email. Any simulated data provided for code testing will contain no UK Biobank participant-level data, real participant identifiers, or real participant information. Researchers who wish to reproduce the full analyses should apply for access through the UK Biobank Access Management System.

Depending on the script used, expected outputs may include processed phenotype files, regression result tables, survival analysis summaries, omics association results, mediation analysis outputs, genetic analysis results, and publication-style figures.

The scripts are provided to document the analytical workflow and support reproducibility among researchers with appropriate data access.

## Data Availability

Individual-level UK Biobank data are not included in this repository. Access to UK Biobank data is available to approved researchers through the official UK Biobank application process.

Derived summary results and example outputs may be shared where permitted by data access policies.

## Code Availability

This repository contains the analytical scripts used for the study. The code will be made publicly available upon publication, and a stable version may be archived in Zenodo or another permanent repository.

## License

This code is released under the MIT License.

See the `LICENSE` file for details.

## Citation

The final citation should be updated after publication.

## Contact

For questions regarding the code or analytical workflow, please contact:

Dr. Zhiqian Huang

Eye & ENT Hospital, and Institute of Science and Technology for Brain-Inspired Intelligence, Fudan University, Shanghai 200031 and 200433, China

[drhuangzhiqian@163.com]


