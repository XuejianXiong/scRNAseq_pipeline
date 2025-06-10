# 🔬 Perturb-seq Analysis Pipeline

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
- pandas, seaborn, matplotlib

## 📈 Key Results

- UMAP plots to visualize cell states by perturbation
- Identification of differentially expressed (DE) genes
- Functional enrichment of DE genes (GO / KEGG pathways)

## ⚙️ Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate perturbseq
```

## 🚀 How to Run

Execute the pipeline step-by-step using the Python scripts under the `scripts/` directory.

```bash
python scripts/01_download_data.py             # Download and extract GSE149383
python scripts/02_preprocessing.py             # Load, filter, and normalize data
python scripts/03_qc.py                        # Perform quality control
python scripts/04_dimred_clustering.py         # Dimensionality reduction + clustering
python scripts/05_perturbation_effects.py      # Differential expression analysis
python scripts/06_biological_analysis.py       # Pathway enrichment analysis
```

Each script reads from `data/` and saves outputs to `figures/` or back into `data/`.

## 📂 Folder Structure

```
perturbseq-pipeline/
├── data/                  # Raw and processed datasets
├── scripts/               # Python scripts for each pipeline step
├── figures/               # Output visualizations
├── environment.yml        # Conda environment definition
├── README.md              # This file
├── .gitignore             # Ignored files/folders
```

## 📘 License

MIT License – feel free to use, adapt, and share.
