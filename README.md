# 🚀 Overview of Pipeline Steps

## 🔹 1. Download Public Perturb-seq Dataset

Goal: Automate fetching and extracting real Perturb-seq data

Tool: wget, tar, gunzip

## 🔹 2. Preprocess scRNA-seq Data + Guide Assignments

Goal: Load 10x-style count matrix and merge with sgRNA annotations

Tool: scanpy, pandas, anndata

## 🔹 3. Quality Control (QC)

Goal: Filter low-quality cells, genes, and empty droplets

Tool: scanpy, thresholding, violin plots

## 🔹 4. Normalization + Dimensionality Reduction

Goal: Normalize, log-transform, run PCA + UMAP

Tool: scanpy.pp functions

## 🔹 5. Clustering and Marker Gene Analysis

Goal: Identify clusters and marker genes

Tool: scanpy.tl.leiden, scanpy.tl.rank_genes_groups

## 🔹 6. Guide-Level Perturbation Effects

Goal: Compare expression profiles of cells with different guides

Tool: Differential expression analysis, volcano plots

## 🔹 7. Biological Interpretation

Goal: Infer pathways affected by perturbations (e.g., GSEA or GO)

Tool: gseapy, scanpy.tl.dendrogram, violin plots

## 🔹 8. Generate Summary Report & Push to GitHub

Goal: Markdown README, summary figures, clean code repo

Tool: matplotlib, seaborn, nbconvert, Git

## 🧪 Example Repo Structure
          perturbseq-pipeline/
          ├── data/                        # Downloaded and processed data
          ├── notebooks/                  # Jupyter notebooks for each step
          │   ├── 01_download_data.ipynb
          │   ├── 02_preprocessing.ipynb
          │   ├── 03_qc.ipynb
          │   ├── 04_dimred_clustering.ipynb
          │   ├── 05_perturbation_effects.ipynb
          │   └── 06_biological_analysis.ipynb
          ├── scripts/                    # Optional .py scripts
          ├── figures/                    # UMAPs, DE results, volcano plots
          ├── README.md                   # Project summary
          ├── environment.yml             # Conda environment
          └── .gitignore
