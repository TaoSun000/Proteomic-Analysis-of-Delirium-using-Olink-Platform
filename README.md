# Proteomic Analysis of Delirium Post-Cardiac Surgery using Olink Platform

## Overview
This repository contains the R code and analysis workflow used in the study investigating protein alterations in patients experiencing delirium following cardiac surgery. The study employs the Olink proteomics platform, involving advanced statistical analysis and visualization techniques to identify significant protein biomarkers and pathways.

## Contents
- **paired_ttest_all.csv**: Results of paired t-tests for all proteins.
- **paired_ttest_significant.csv**: List of proteins significantly altered at 5% FDR.
- **top20_volcanoplot.pdf**: Volcano plot visualization of top 20 significant proteins.
- **top5_protein_boxplot.pdf**: Boxplot visualization for the top 5 proteins.
- **pca_analysis.pdf**: PCA analysis plot distinguishing case and control groups.
- **GSEA_results.csv**: Results from Gene Set Enrichment Analysis (GSEA).
- **GSEA_heatmap.pdf**: Heatmap visualization of significant pathways identified by GSEA.
- **ORA_results.csv**: Results from Over-Representation Analysis (ORA).
- **ORA_heatmap.pdf**: Heatmap visualization of significant pathways identified by ORA.
- **ROC_curve_FKBP1B.pdf**: ROC curve for the top protein (FKBP1B).

## Data Files (replace paths accordingly)
- Batch1 and Batch2 raw data files in NPX format.
- Updated status information for matched case-control pairs.

## Dependencies
R packages required:
- `OlinkAnalyze`
- `dplyr`
- `tidyr`
- `stringr`
- `ggplot2`
- `clusterProfiler`
- `org.Hs.eg.db`
- `pheatmap`
- `ggrepel`
- `survival`
- `pROC`

Install missing dependencies using:
```R
install.packages(c("OlinkAnalyze", "dplyr", "tidyr", "stringr", "ggplot2", 
                   "clusterProfiler", "org.Hs.eg.db", "pheatmap", "ggrepel", 
                   "survival", "pROC"))
```

## Running the Analysis
Clone this repository and execute the R script (`proteomics_analysis.R`) in an environment configured with the necessary packages listed above.

## Contact
For any questions or feedback, please contact:
- Tao Sun, MD, PhD, MSc  
- Email: taosun618@gmail.com

---

© 2025 Tao Sun. All rights reserved.
