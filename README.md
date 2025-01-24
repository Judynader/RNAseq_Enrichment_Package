# RNAseq_Enrichment_Package

This package is designed for RNA-seq data analysis, DEG identification, and enrichment visualization.

## Collaboration
This project was developed as part of a collaborative effort:
- **Yael Campoy Fernández**: Main developer responsible for writing the R code.
- **Judy Alshedah**: GitHub management, version control, and documentation.

## Features
- Import RNA-seq count data and filter low-expressed genes.
- Identify differentially expressed genes (DEGs) using edgeR.
- Perform enrichment analysis (GO and KEGG pathways) using clusterProfiler.
- Export results to Excel and visualize enrichment results (dotplot, cnetplot, treeplot).
## Installation
You can install the package from GitHub using devtools:
```r
install.packages("devtools")
devtools::install_github("JudyNader/RNAseq_Enrichment_Package")
