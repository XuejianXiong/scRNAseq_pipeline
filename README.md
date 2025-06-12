# 🔬 Single-cell RNA-seq Analysis Pipeline

This project analyzes public CROP-seq data from A549 cells (GSE149383) to study transcriptional effects of CRISPR-based gene perturbations using a reproducible single-cell RNA-seq pipeline built in Python.

## 📊 Dataset

**Citation:**  
Replogle et al. (2020). *Direct capture of CRISPR guides enables scalable, multiplexed, and multi-omic Perturb-seq*. Cell.  
GEO Accession: [GSE149383](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149383)

- Cell line: A549 (lung cancer)
- Platform: CRISPRi + 10x Genomics
- Perturbations: Targeted gene knockdowns via CROP-seq

## 🛠️ Tech Stack

- Python 3.10
- [Scanpy](https://scanpy.readthedocs.io/)
- [gseapy](https://gseapy.readthedocs.io/)
- pandas, seaborn, matplotlib, numpy, anndata
- python-igraph, leidenalg

## 📈 Key Results

- UMAP plots to visualize cell states by perturbation
- Identification of differentially expressed (DE) genes
- Functional enrichment of DE genes (GO / KEGG pathways)

## ⚙️ Installation

Create the python virtural environment and install required python packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## 🚀 How to Run

Execute the pipeline step-by-step using the Python scripts under the `scripts/` directory.

```bash
python scripts/00_setup.sh                     # Setup directories
python scripts/01_download_data.sh             # Download and extract GSE149383
python scripts/02_preprocessing.py             # Load, filter, and normalize data
python scripts/03_qc.py                        # Perform quality control
python scripts/04_normalization_dimred.py      # Normalization + dimensionality reduction
python scripts/05_clustering.py                # Clustering
python scripts/06_DE.py                        # Differential expression analysis
python scripts/07_GSEA.py                      # Pathway enrichment analysis
```

Each script reads from `data/` and saves outputs to `results/` and `figures/`.

## 📂 Folder Structure

```
perturbseq-pipeline/
├── data/                  # Raw and processed datasets
├── scripts/               # Python scripts for each pipeline step
├── figures/               # Output visualizations
├── results/               # Output data files
├── requirements.txt       # Python packages
├── README.md              # This file
├── .gitignore             # Ignored files/folders
```

## 📘 License

MIT License – feel free to use, adapt, and share.
