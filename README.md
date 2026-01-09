# REHAB-HF Proteomics Analysis

Comprehensive proteomics analysis of protein biomarkers and rehabilitation outcomes in heart failure.

**Study:** Secondary Analysis of REHAB-HF (Rehabilitation and Exercise Training After Hospitalization [NCT01508650])

**Primary Outcomes:** 
- SPPB (Short Physical Performance Battery)
- 6MWD (Six-Minute Walk Distance)

## Overview

This repository contains the complete analysis pipeline for identifying and validating protein biomarkers that modify the effectiveness of rehabilitation intervention on physical performance in heart failure patients. The analysis employs multiple statistical and machine learning approaches to screen 92 Olink proteomics biomarkers and identify those most predictive of differential treatment response.

### Protein Biomarkers

- **92 biomarkers**
- **5 key proteins** selected for primary analysis using individual regression analyses, LASSO screening with strong hierarchy, and Bayesian additive regression tree interaction analysis :
  - ALCAM
  - Gal-4
  - PECAM-1
  - GP6
  - ST2

## Repository Structure

RehabHFproteomics/
├── README.md                          # This file
├── REHAB-HF_proteomicsanalysis010726_GitHub.R 
├── Results/ #contains all tables, figures, and additional output
│   ├── RegressionAnalyses/ #contains summary output of individual regression analyses

## Getting Started

### Requirements

- R ≥ 4.3.0
- RStudio (recommended)
- Quarto or R Markdown

### Required R Packages

```r
library(tidyverse)
library(XLConnect)
library(readxl)
library(janitor)
library(gt)
library(gtsummary)
library(reshape2)
library(sjPlot)
library(ggforce)
library(gridExtra)
library(glinternet)
library(bartMachine)
library(optmatch)
library(rpart)
library(rpart.plot)
library(lsmeans)
library(treeClust)
library(pdfCluster)
library(mediation)
library(dplyr)
library(ggcorrplot)
library(HDtest)
library(WGCNA)
library(RColorBrewer)
library(ggrepel)
library(forcats)
library(randomForest)
```

## Data Dictionary

### Demographic/Clinical Variables

| Variable | Description | Values |
|----------|-------------|--------|
| `age` | Age in years | continuous |
| `sex` | Biological sex | 1=Female, 0=Male |
| `race___4` | Race/ethnicity | 1=White, 0=Non-White |
| `hf_cat` | HF classification | 1=HFpEF (EF≥45%), 0=HFrEF (EF<45%) |
| `creatinine_value` | Serum creatinine | mg/dL |
| `egfr` | Estimated GFR | mL/min/1.73m² (2021 CKD-EPI) |

### Physical Performance Outcomes

| Variable | Description | Range |
|----------|-------------|-------|
| `bl_sppb` | Baseline SPPB score | 0-12 (higher=better) |
| `fu_sppb` | 12-week follow-up SPPB | 0-12 |
| `bl_smw` | Baseline 6MWD | meters |
| `fu_smw` | 12-week follow-up 6MWD | meters |

### Intervention

| Variable | Description | Values |
|----------|-------------|--------|
| `intervention_1_control_0` | Treatment assignment | 1=Rehabilitation, 0=Attention Control |

## Analysis Sections

### Section 1: Baseline Characteristics (Table 1)
Summary statistics stratified by intervention group with baseline tests for group differences.

### Section 2: Linear Regression (Table 2)
Primary analysis testing associations between protein expression and 12-week outcomes.
- Unadjusted and adjusted models
- Tests for treatment × protein interaction

### Sections 3-4: Figures 1-3
Visualizations of protein correlations and protein-outcome relationships.

### Section 5: Propensity Score Matching & Matching Trees
Causal inference using 1:1 matching on propensity scores and protein values. Matching trees identify subgroups with different treatment effects.

### Section 6: WGCNA
Identifies co-expression modules of proteins and their associations with clinical outcomes.

### Section 9: Individual Effect Modification Screening
Screens all 92 proteins for treatment interactions (unadjusted and FDR-corrected p-values).

### Sections 10-11: LASSO & BART Screening
Machine learning variable selection:
- **LASSO** (glinternet): L1 regularization with hierarchical constraints
- **BART**: Tree-based ensemble method with interaction importance

### Section 12: Mediation Analysis
Tests whether changes in protein expression mediate intervention effects on outcomes.

## Output Files

### Main Results
- `Table1_Baselinecharacteristics.html` - Demographic summary
- `Table2_sppb.csv` - SPPB regression results
- `Table2_6mwd.csv` - 6MWD regression results

### Figures
- `Figure1_SpearmanCorrelation.pdf` - Protein correlation heatmap
- `Figure2_ChangeinSPPB.pdf` - Protein × intervention on SPPB change
- `Figure3_Changein6MWD.pdf` - Protein × intervention on 6MWD change
- `Figure4a_Matchingtree_SPPB.pdf` - Matching tree (SPPB)
- `Figure4b_Matchingtree_6MWD.pdf` - Matching tree (6MWD)

### Supplemental Files
- `SuppFigure2_InteractionImportance_sppb.pdf` - BART importance (SPPB)
- `SuppFigure2_InteractionImportance_6mwd.pdf` - BART importance (6MWD)
- `SuppFigure3_CVPlot_*.pdf` - Cross-validation curves
- `SuppFigure4_Proteinchange.pdf` - Protein changes by intervention
- `protein_modules.csv` - WGCNA module assignments

## Contact

Questions or issues? Please open an issue on GitHub or contact the corresponding author.

**Repository:** https://github.com/sbruce23/RehabHFproteomics

This README provides comprehensive documentation for your GitHub repository with clear structure, installation instructions, and complete analysis overview.
