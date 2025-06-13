# 🔬 Single-cell RNA-seq Analysis Pipeline (Python + R)

This project presents a complete and reproducible **single-cell RNA-seq (scRNA-seq) analysis pipeline** implemented in both **Python** and **R**. It uses public **CROP-seq data** (CRISPRi + 10x Genomics) from A549 lung cancer cells to explore transcriptional effects of gene perturbations at the single-cell level.

---

## 📊 Dataset

**Study:**  
Replogle et al. (2020). *Direct capture of CRISPR guides enables scalable, multiplexed, and multi-omic Perturb-seq*. **Cell**  
**GEO Accession:** [GSE149383](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149383)

- Cell line: A549 (lung adenocarcinoma)
- Technology: CRISPRi + 10x Genomics
- Platform: CROP-seq
- Objective: Identify transcriptional changes in response to gene knockdowns

---

## 🧰 Tech Stack

### Python version
- Python 3.13.3
- [Scanpy](https://scanpy.readthedocs.io/)
- [gseapy](https://gseapy.readthedocs.io/)
- pandas, seaborn, matplotlib, numpy, anndata
- python-igraph, leidenalg

### R version
- R 4.5.0
- [Seurat](https://satijalab.org/seurat/)
- dplyr, ggplot2, patchwork
- fgsea, clusterProfiler, org.Hs.eg.db
- Bioconductor & tidyverse ecosystem

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
python python_scripts/00_setup.sh                     # Setup directories
python python_scripts/01_download_data.sh             # Download and extract GSE149383
python python_scripts/02_preprocessing.py             # Load, filter, and merge datasets
python python_scripts/03_qc.py                        # Perform quality control
python python_scripts/04_normalization_dimred.py      # Normalize and run PCA/UMAP
python python_scripts/05_clustering.py                # Clustering (Leiden)
python python_scripts/06_DE.py                        # Differential expression
python python_scripts/07_GSEA.py                      # Pathway enrichment (GO/KEGG)
```

### 🟣 R Pipeline

```bash
source("install_packages.R")
```

Run R scripts in RStudio or VS Code:

```bash
source("R_scripts/00_setup.R")                  # Set up directories
source("R_scripts/01_download_data.R")          # Download and extract GSE149383
source("R_scripts/02_preprocessing.R")          # Merge datasets with metadata
source("R_scripts/03_qc.R")                     # Perform quality control
source("R_scripts/04_normalization_dimred.R")   # Normalize and run PCA/UMAP
source("R_scripts/05_clustering.R")             # Clustering
source("R_scripts/06_DE.R")                     # DE analysis using Seurat
source("R_scripts/07_GSEA.R")                   # Enrichment analysis using fgsea
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
└── data/
|   └── GSE149383/         # Raw and processed datasets
├── python_scripts/        # Python scripts for each pipeline step
├── R_scripts/             # R scripts for each pipeline step

```

---

## 🧪 Key Results

UMAP visualization of cell states under gene perturbations

Identification of differentially expressed genes (DEGs)

Functional enrichment (GO/KEGG) of DEGs

Side-by-side implementation in Python and R for reproducibility and flexibility


---

## 📘 License

MIT License – feel free to use, adapt, and share.
