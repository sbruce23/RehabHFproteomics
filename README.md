# REHAB-HF Proteomics Analysis

Comprehensive proteomics analysis of protein biomarkers and rehabilitation outcomes in heart failure.

Full analysis available at https://sbruce23.github.io/RehabHFproteomics/REHAB-HF_ProteomicsMarkdown_cleaned_011326.html

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

## Contact

Questions or issues? Please open an issue on GitHub or contact the corresponding author.

**Repository:** https://github.com/sbruce23/RehabHFproteomics
