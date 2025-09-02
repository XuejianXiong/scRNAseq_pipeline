# 🔬 Single-cell RNA-seq Analysis Pipeline  

[![Version](https://img.shields.io/badge/version-2.0-blue.svg)](https://github.com/yourname/scRNAseq_pipeline/releases/tag/v2.0)  
[![Python](https://img.shields.io/badge/python-3.13+-brightgreen.svg)](https://www.python.org/)  
[![R](https://img.shields.io/badge/R-4.5.0+-purple.svg)](https://www.r-project.org/)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  

## 📌 Version History

- **[v2.0](https://github.com/yourname/scRNAseq_pipeline/tree/main)** – Modular Python pipeline + legacy R (**current main branch**)  
- **[v1.0](https://github.com/yourname/scRNAseq_pipeline/tree/feature_qc)** – Original Python + R pipeline (non-modular)  

---

## **Version 2 – Modular Python + Legacy R**

This project provides a complete and reproducible **single-cell RNA-seq (scRNA-seq) analysis pipeline** implemented in both **Python** and **R**.  

- **Version 1 (feature_qc branch):** Original pipeline in **Python and R**, without modular design.  
- **Version 2 (main branch):** **Modular, robust Python pipeline** with **R scripts unchanged from Version 1**.  

### **Datasets processed in this pipeline**
- **CROP-seq data** (CRISPRi + 10x Genomics) from A549 lung cancer cells – *GSE149383*  
- **Retina datasets**:  
  - *SRA559821* (from [PanglaoDB](https://panglaodb.se/))  
  - *GSE137537* – from *"Single-cell Transcriptomic Atlas of the Human Retina Identifies Cell Types Associated with Age-Related Macular Degeneration"*  

---

## 📊 Dataset

### 1. CROP-seq A549 Perturbation
- **Study:**  Replogle et al. (2020). *Direct capture of CRISPR guides enables scalable, multiplexed, and multi-omic Perturb-seq*. **Cell**  
- **GEO Accession:** [GSE149383](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149383)
- Cell line: A549 (lung adenocarcinoma)
- Technology: CRISPRi + 10x Genomics
- Platform: CROP-seq
- Objective: Identify transcriptional changes in response to gene knockdowns

### 2. Retina scRNA-seq Datasets
- **SRA559821 (PanglaoDB)** – Reference retina dataset for cell type annotation  
- **GSE137537** – *Human Retina Transcriptomic Atlas (Age-related Macular Degeneration)*  
- Objective: Identify and compare retina cell populations and disease-associated transcriptional signatures  

---

## 🧰 Tech Stack

### **Python (Version 2 – Modular)**
- Python 3.13.3  
- [Scanpy](https://scanpy.readthedocs.io/) for scRNA-seq analysis  
- [gseapy](https://gseapy.readthedocs.io/) for pathway enrichment  
- pandas, numpy, matplotlib, seaborn, anndata  
- python-igraph, leidenalg  

### **R (Unchanged from Version 1)**
- R 4.5.0  
- [Seurat](https://satijalab.org/seurat/), SeuratObject  
- dplyr, ggplot2, patchwork, readr, tibble, Matrix  
- fgsea, msigdbr, pheatmap, knitr  

> **Note:** R scripts remain from Version 1 and are fully functional, but **not modularized yet**. Future updates will align the R workflow with the robust Python structure.

---

## 🚀 How to Run the Pipelines

> Each version can be run independently. Output folders and filenames are standardized.

### 🔷 Python Pipeline

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Then run each step:

```bash
./python_scripts/01_download_data_cropseq.sh                                              # Download CROP-seq data
python python_scripts/01_download_GEOretina.py                                            # Download retina GSE137537 data
python python_scripts/01_convert_panglao_to_10x.py                                        # Convert Panglaodb data to 10x format
python python_scripts/02_preprocessing.py cropseq python_scripts/config.yaml              # Load, filter, and merge datasets
python python_scripts/03_qc.py cropseq python_scripts/config.yaml                         # Perform quality control
python python_scripts/04_normalization_dimred.py cropseq python_scripts/config.yaml       # Normalize and run PCA/UMAP
python python_scripts/05_clustering.py cropseq python_scripts/config.yaml                 # Clustering (Leiden)
python python_scripts/06_DE.py cropseq python_scripts/config.yaml                         # Differential expression
python python_scripts/07_GSEA.py cropseq python_scripts/config.yaml                       # Pathway enrichment (GO/KEGG)
```

### 🟣 R Pipeline

```bash
source("install_packages.R")
```

Run R scripts in RStudio or VS Code:

```bash
./R_scripts/00_setup.sh                           # Set up directories
./R_scripts/01_download_data.sh                   # Download data
source("R_scripts/02_preprocessing.R")            # Merge datasets with metadata
source("R_scripts/03_qc.R")                       # Perform quality control
source("R_scripts/04_normalization_dimred.R")     # Normalize and run PCA/UMAP
source("R_scripts/05_clustering.R")               # Clustering
source("R_scripts/06_DE.R")                       # DE analysis using Seurat
source("R_scripts/07_GSEA.R")                     # Enrichment analysis using fgsea
```

---

## 📂 Folder Structure

```
scRNAseq_pipeline/
├── README.md              # This file
├── .gitignore             # Ignored files/folders
├── requirements.txt       # Python packages
├── install_packages.R     # R packages
|
├── figures/               # Output visualizations
├── results/               # Output data files
├── data/                  # Input data files
├── R_scripts/             # R scripts for each pipeline step
├── python_scripts/        # Python scripts for each pipeline step

```

---

## 🧪 Key Results

UMAP visualization of perturbation and retina cell states

Identification of differentially expressed genes (DEGs) across multiple datasets

Functional enrichment (GO/KEGG) of DEGs

Modular, maintainable design in Python (Version 2)

Legacy R scripts kept for reproducibility (Version 1)

---

## 📘 License

MIT License – feel free to use, adapt, and share.
