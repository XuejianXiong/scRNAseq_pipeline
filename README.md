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
- Python 3.10
- [Scanpy](https://scanpy.readthedocs.io/)
- [gseapy](https://gseapy.readthedocs.io/)
- pandas, seaborn, matplotlib, numpy, anndata
- python-igraph, leidenalg

### R version
- R (≥ 4.2.0)
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
python scripts/00_setup.sh                     # Setup directories
python scripts/01_download_data.sh             # Download and extract GSE149383
python scripts/02_preprocessing.py             # Load, filter, and merge datasets
python scripts/03_qc.py                        # Perform quality control
python scripts/04_normalization_dimred.py      # Normalize and run PCA/UMAP
python scripts/05_clustering.py                # Clustering (Leiden)
python scripts/06_DE.py                        # Differential expression
python scripts/07_GSEA.py                      # Pathway enrichment (GO/KEGG)
```

### 🟣 R Pipeline

Open R/sc_pipeline.Rproj or run scripts manually in RStudio or VS Code:

```bash
source("R/00_setup.R")               # Set up directories
source("R/01_load_data.R")           # Load 10x data from nested folders
source("R/02_preprocessing.R")       # Merge datasets with metadata
source("R/03_qc.R")                  # Quality control and filtering
source("R/04_dimred_clustering.R")   # PCA, UMAP, clustering
source("R/05_DE_analysis.R")         # DE analysis using Seurat
source("R/06_GSEA.R")                # Enrichment analysis using fgsea
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
    └── GSE149383/         # Raw and processed datasets
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
