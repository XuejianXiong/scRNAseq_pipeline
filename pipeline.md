# 🚀 Steps of the Pertub-seq Analysis Pipeline

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
          ├── scripts/                     # .py scripts for each step
          │   ├── 01_download_data.py
          │   ├── 02_preprocessing.py
          │   ├── 03_qc.py
          │   ├── 04_dimred_clustering.py
          │   ├── 05_perturbation_effects.py
          │   └── 06_biological_analysis.py
          ├── figures/                    # UMAPs, DE results, volcano plots
          ├── results/                    # processed data
          ├── README.md                   # Project summary
          ├── environment.yml             # Conda environment
          └── .gitignore



Biological Interpretation of QC and HVG Plots
Quality Control Metrics
Number of Genes Detected per Cell:
This metric indicates the complexity of the transcriptome captured in each cell. Cells with very low gene counts often represent low-quality or dying cells and should be filtered out to avoid noise.

Total Counts per Cell:
The total number of transcripts detected per cell can reflect sequencing depth or capture efficiency. Cells with abnormally high counts might be doublets or multiplets.

Percent Mitochondrial Gene Counts:
High mitochondrial RNA percentage often signals stressed or dying cells, as mitochondrial transcripts tend to increase when cells are undergoing apoptosis or other damage.

Filtering based on these metrics ensures downstream analysis focuses on biologically meaningful cells.

Highly Variable Genes (HVGs)
HVGs are genes whose expression varies significantly across cells beyond technical noise. These genes often represent key biological signals such as cell identity, states, or responses to stimuli.

Focusing on HVGs during downstream analyses (like clustering or dimensionality reduction) enhances the detection of meaningful cell subpopulations and biological processes.

Biologically, HVGs might include transcription factors, signaling molecules, or other genes driving cellular heterogeneity.

Dimensionality Reduction (PCA, UMAP)
PCA reduces data complexity while retaining major sources of variation, helping identify principal axes of biological variability.

UMAP further embeds cells into a low-dimensional space that preserves local and global relationships, allowing intuitive visualization of cell clusters or trajectories.

Together, these techniques reveal the structure of the cellular landscape, highlight subpopulations, and suggest biological distinctions such as cell types or states.





# GSE149383 Single-cell RNA-seq Analysis Pipeline (10x Genomics)

## 0. SETUP DIRECTORIES
# shell: 00_setup.sh

mkdir -p data/GSE149383/raw
mkdir -p data/GSE149383/extracted
mkdir -p results

## 1. DOWNLOAD DATA
# shell: 01_download.sh

cd data/GSE149383/raw

# Download raw tar archive (alternatively use wget or curl)
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149383/suppl/GSE149383_RAW.tar

tar -xvf GSE149383_RAW.tar -C ../extracted
cd ../extracted

# Extract relevant samples (example: PC9D0 and PC9D3Erl)
for f in GSM3972651_PC9D0_filtered_feature_bc_matrices.tar.gz GSM3972652_PC9D3Erl_filtered_feature_bc_matrices.tar.gz; do
  tar -xvzf $f
  rm $f
done

## 2. LOAD AND MERGE DATA
# Python: 02_preprocessing.py

import scanpy as sc
import anndata
import os
import pandas as pd

# Load filtered matrices from both samples
data_dir = "data/GSE149383/extracted"
samples = {
    "PC9D0": os.path.join(data_dir, "filtered_feature_bc_matrices", "hg19"),
    "PC9D3Erl": os.path.join(data_dir, "filtered_feature_bc_matrices_2", "hg19")
}

adatas = []
for name, path in samples.items():
    adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    adata.obs['sample'] = name
    adatas.append(adata)

# Merge datasets
adata_merged = adatas[0].concatenate(adatas[1:], batch_key="sample")

# Save merged AnnData
adata_merged.write("results/merged_data.h5ad")

## 3. PREPROCESSING
# Python: 03_preprocess.py

adata = sc.read("results/merged_data.h5ad")

# Basic QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Save filtered object
adata.write("results/filtered_data.h5ad")

## 4. DIMENSIONALITY REDUCTION AND CLUSTERING
# Python: 04_analysis.py

adata = sc.read("results/filtered_data.h5ad")

# PCA + neighbors + clustering
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Save clustering result
adata.write("results/clustered_data.h5ad")

## 5. VISUALIZATION
# Python: 05_plot.py

import seaborn as sns
import matplotlib.pyplot as plt

adata = sc.read("results/clustered_data.h5ad")

# UMAP by sample
sc.pl.umap(adata, color=["sample"], save="_sample.png")

# UMAP by cluster
sc.pl.umap(adata, color=["leiden"], save="_clusters.png")

# Top marker genes per cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, save="_markers.png")

## 6. BIOLOGICAL INTERPRETATION
# Biological idea: PC9D0 (baseline) vs PC9D3Erl (3 days Erlotinib)
# Expect to see treatment-specific gene expression shifts.
# Marker genes will help identify pathway shifts (e.g., EGFR signaling).

# Enrichment analysis (optional)
# Use gseapy or other tools on marker gene lists

# END
