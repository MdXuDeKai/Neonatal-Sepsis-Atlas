# Neonatal-Sepsis-Atlas

**Single-Cell Immunometabolic Profiling Links High Lactate-Metabolizing Monocytes to T Cell Stress Programs in Preterm Neonatal Sepsis**

## 📌 Overview
This repository is dedicated to the analysis code and scripts for our study investigating the immunometabolic dysregulation in preterm infants with late-onset sepsis (LOS). 

By leveraging paired single-cell RNA sequencing (scRNA-seq) and bulk transcriptomic validation, we identified a novel monocyte subpopulation—**High Lactate-Metabolizing Inflammatory Monocytes (HLMs)**—and mapped their putative regulatory crosstalk with T cell stress programs.

## 🚧 Code Availability
**The complete analysis pipeline, processing scripts, and figure reproduction code are currently being finalized and organized.** All code will be made publicly available in this repository immediately upon the **acceptance/publication** of the manuscript.

### 📢 For Reviewers and Editors
If you are a reviewer or editor requiring early access to the code for assessment purposes, or if you have specific questions regarding the methodology, please do not hesitate to contact the corresponding author directly. We are happy to share the scripts upon reasonable request.

📧 **Contact:** [Insert Your Email Here]

## 📂 Data Availability
The datasets analyzed in this study are publicly available in the Gene Expression Omnibus (GEO):
* **scRNA-seq Discovery Cohort:** [GSE236099](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236099)
* **Bulk Validation Cohort:** [GSE25504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25504)

## 🛠️ Key Analysis Tools
The analysis was primarily performed using **R (v4.5.1)**. Key packages include:
* **Seurat** (v4.4.0) - Quality control, clustering, and dimensionality reduction.
* **Harmony** (v1.2.0) - Batch effect correction.
* **MAST** (v1.22.0) - Differential expression analysis.
* **CellChat** (v1.6.1) & **NicheNet** (v1.1.0) - Ligand-receptor communication inference.
* **clusterProfiler** (v4.4.4) - Functional enrichment analysis.

---
*Copyright © 2026. All rights reserved.*
