# -----------------------------
# Step 4: Normalization and Dimensionality Reduction
# -----------------------------

import scanpy as sc
import numpy as np
import pandas as pd
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# User-Adjustable Parameters
# -----------------------------
RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("04_log.txt")
INPUT_FILE = RESULTS_DIR / "03_filtered_data.h5ad"
OUTPUT_FILE = RESULTS_DIR / "04_norm_dr_data.h5ad"
HTML_REPORT = RESULTS_DIR / "04_report.html"

FIG_NAME1 = "04_hvg_plot.png"
FIG_FILE2 = FIGURE_DIR / "04_hvg_colored.png"
FIG_NAME3 = "04_pca_treatment.png"
FIG_NAME4 = "04_umap_treatment.png"

N_TOP_HVG = 2000
SCALE_MAX = 10
N_PCS = 40
N_NEIGHBORS = 15
PALETTE_NAME = "Set2"
PALETTE_SIZE = 10
# -----------------------------


logging.info("Step 4 started: Normalization and Dimensionality Reduction")

# -----------------------------
# Define unified color palette
# -----------------------------
treatment_palette = sns.color_palette(PALETTE_NAME, n_colors=PALETTE_SIZE)
treatment_categories = None

# -----------------------------
# Helper Function: Custom HVG Plot
# -----------------------------
def plot_highly_variable_genes(adata, save_path):
    if 'highly_variable' not in adata.var.columns:
        raise ValueError("HVGs not computed. Run `sc.pp.highly_variable_genes` first.")
    
    means = adata.var['means']
    dispersions = adata.var['dispersions']
    hvg_mask = adata.var['highly_variable']

    plt.figure(figsize=(6, 4))
    plt.scatter(means[~hvg_mask], dispersions[~hvg_mask], c='lightgray', s=10, label='Not HVG', alpha=0.5)
    plt.scatter(means[hvg_mask], dispersions[hvg_mask], c='red', s=10, label='HVG', alpha=0.8)
    plt.xlabel('Mean expression')
    plt.ylabel('Dispersion')
    plt.title('Highly Variable Genes')
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

# -----------------------------
# Load Data
# -----------------------------
adata = sc.read(INPUT_FILE)
logging.info(f"Loaded filtered data: {adata.n_obs} cells, {adata.n_vars} genes")

# Ensure 'treatment' is categorical and assign consistent colors
adata.obs["treatment"] = adata.obs["treatment"].astype("category")
treatment_categories = adata.obs["treatment"].cat.categories
treatment_colors = dict(zip(treatment_categories, treatment_palette[:len(treatment_categories)]))
adata.uns["treatment_colors"] = [treatment_colors[cat] for cat in treatment_categories]
logging.info(f"Set treatment categories and colors: {treatment_colors}")

# -----------------------------
# Normalization and HVG Detection
# -----------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
logging.info("Normalized total counts per cell to 10,000")

sc.pp.log1p(adata)
adata.raw = adata 
logging.info("Applied log1p transformation")

sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_HVG, flavor="seurat", subset=True)
n_hvg = np.sum(adata.var["highly_variable"])
logging.info(f"Selected {n_hvg} highly variable genes")

sc.pl.highly_variable_genes(adata, save=f"_{FIG_NAME1}", show=False)
plot_highly_variable_genes(adata, save_path=FIG_FILE2)
logging.info("Saved HVG plots")

# -----------------------------
# Scaling, PCA, and UMAP
# -----------------------------
sc.pp.scale(adata, max_value=SCALE_MAX)
logging.info(f"Scaled the data to unit variance with max value {SCALE_MAX}")

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color="treatment", palette=treatment_colors, save=f"_{FIG_NAME3}", show=False)
logging.info("Performed PCA and saved plot colored by treatment")

sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
sc.tl.umap(adata)
sc.pl.umap(adata, color="treatment", palette=treatment_colors, save=f"_{FIG_NAME4}", show=False)
logging.info(f"Computed neighbors (n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}) and UMAP, saved plot colored by treatment")

# -----------------------------
# Save Processed Data
# -----------------------------
adata.write(OUTPUT_FILE)
logging.info("Saved normalized and reduced AnnData object")

# -----------------------------
# Generate HTML Report
# -----------------------------
figures = sorted(Path(FIGURE_DIR).glob("*04_*.png"))

with open(HTML_REPORT, "w") as f:
    f.write("<html><head><title>Step 4: Normalization & Dimensionality Reduction</title></head><body>\n")
    f.write("<h1>Step 4 Report: Normalization, HVG, PCA, UMAP</h1>\n")
    for fig in figures:
        f.write(f"<h3>{fig.name}</h3>\n")
        f.write(f'<img src="../{fig}" width="700"><br><br>\n')
    f.write("</body></html>")
logging.info(f"Generated HTML report at {HTML_REPORT}")


print("✅ Step 4 complete: Normalization, HVG selection, PCA & UMAP done and saved.")
