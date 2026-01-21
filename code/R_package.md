## R package installation guide (this project)

This document lists the R packages used across this repository and how to install them (CRAN / Bioconductor / GitHub). It is intended to help reviewers/readers reproduce the analysis environment.

> **Recommended**: run the analysis in **RStudio** on Windows for more robust handling of fonts/encoding.

---

## 1) Prerequisites

- **R**: >= 4.0.0 (recommended R >= 4.3)
- **Rtools (Windows)**: recommended if you need to compile packages from source (GitHub installs may trigger compilation).
- **Internet access**: some packages are installed from CRAN/Bioconductor/GitHub.

---

## 2) One-shot installation (recommended)

Copy-paste the following into an R session started at the project root.

```r
## ------------------------------------------------------------
## Install all packages used in this repository (non-interactive)
## ------------------------------------------------------------

## CRAN
cran_pkgs <- c(
  "Seurat", "hdf5r", "harmony",
  "dplyr", "tidyr", "tibble", "readr", "data.table",
  "ggplot2", "patchwork", "scales", "viridis", "RColorBrewer",
  "ggrepel", "ggpubr",
  "pheatmap", "corrplot", "Hmisc",
  "igraph", "ggraph",
  "remotes", "devtools",
  "NMF",
  "openxlsx"   # optional (only needed for Excel export)
)
install.packages(setdiff(cran_pkgs, rownames(installed.packages())), dependencies = TRUE)

## Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c(
  "MAST", "DESeq2", "edgeR", "limma",
  "clusterProfiler", "fgsea", "msigdbr",
  "org.Hs.eg.db", "AnnotationDbi",
  "ComplexHeatmap", "circlize", "Biobase",
  "GEOquery", "GSVA"
)
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())),
                     ask = FALSE, update = FALSE, dependencies = TRUE)

## GitHub packages used by this project
## - DoubletFinder (doublet detection)
## - presto (optional, accelerates FindAllMarkers)
## - CellChat (cell–cell communication)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", dependencies = TRUE, upgrade = "never")
remotes::install_github("immunogenomics/presto", dependencies = TRUE, upgrade = "never")
remotes::install_github("sqjin/CellChat", dependencies = TRUE, upgrade = "never")

## NicheNet
## Some scripts install nichenetr via BiocManager; if needed, run:
## BiocManager::install("nichenetr", ask = FALSE, update = FALSE)
```

---

## 3) Package list by analysis module

### Core single-cell processing (Seurat objects, QC, integration)

- Seurat
- hdf5r
- dplyr, tidyr, data.table, tibble, readr
- ggplot2, patchwork, scales, viridis, RColorBrewer, ggrepel
- harmony

### Doublet detection / speed-ups

- DoubletFinder (GitHub)
- presto (GitHub; optional, speeds up `FindAllMarkers`)

### Paired differential expression / pseudobulk validation

- MAST (Bioconductor)
- DESeq2 (Bioconductor)
- edgeR, limma (Bioconductor; used for pseudobulk in NicheNet pipeline)

### Enrichment / GSEA

- clusterProfiler (Bioconductor)
- fgsea (Bioconductor)
- msigdbr (Bioconductor/CRAN depending on your setup; installed via BiocManager above for convenience)
- org.Hs.eg.db (Bioconductor)

### NicheNet (ligand–receptor–target inference)

- nichenetr (installed via BiocManager in project scripts)
- ComplexHeatmap (Bioconductor), circlize (Bioconductor)
- igraph, ggraph

### CellChat (communication visualization / mechanistic validation)

- CellChat (GitHub)
- ComplexHeatmap, circlize, pheatmap

### External bulk validation (GSE25504)

- GEOquery (Bioconductor)
- limma (Bioconductor)
- GSVA (Bioconductor)
- AnnotationDbi, org.Hs.eg.db
- ComplexHeatmap, circlize
- corrplot

### Optional exports

- openxlsx (optional; only needed when exporting Excel files)

---

## 4) Sanity check

After installation, you can quickly verify everything loads:

```r
pkgs <- c(
  "Seurat","hdf5r","harmony","dplyr","tidyr","data.table","tibble","readr",
  "ggplot2","patchwork","scales","viridis","RColorBrewer","ggrepel","ggpubr",
  "pheatmap","corrplot","Hmisc","igraph","ggraph",
  "MAST","DESeq2","edgeR","limma",
  "clusterProfiler","fgsea","msigdbr","org.Hs.eg.db","AnnotationDbi",
  "ComplexHeatmap","circlize",
  "GEOquery","GSVA",
  "nichenetr","CellChat","DoubletFinder"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  message("Missing packages:\n- ", paste(missing, collapse = "\n- "))
} else {
  message("All packages are available.")
}
sessionInfo()
```

---

## 5) Notes specific to this repository

- Some scripts include a local `setwd("E:/GBA465/败血症单细胞")`. Please change it to your local project path, or remove `setwd()` and set the working directory in RStudio.
- For CellChat installation details, see `Part3_Step3_Install_CellChat.R` (supports local zip install or GitHub install).
