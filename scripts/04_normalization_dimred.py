import scanpy as sc
import numpy as np
import pandas as pd
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os

# -----------------------------
# Step 4: Normalization and Dimensionality Reduction
# -----------------------------

# -----------------------------
# User-Adjustable Parameters
# -----------------------------
# --- USER PARAM
FIGURE_DIR = "figures"
RESULTS_DIR = "results"
N_TOP_HVG = 2000
SCALE_MAX = 10
N_PCS = 40
N_NEIGHBORS = 15
PALETTE_NAME = "Set2"
PALETTE_SIZE = 10
LOG_FILE = f"{RESULTS_DIR}/04_log.txt"
# --- END USER PARAM

# Setup: logging, figure directory, and result paths
sc.settings.verbosity = 2
sc.settings.autoshow = False
sc.settings.figdir = FIGURE_DIR
Path(FIGURE_DIR).mkdir(parents=True, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format="%(asctime)s %(message)s")
logging.info("Step 4 started: Normalization and Dimensionality Reduction")

# -----------------------------
# Define unified color palette
# -----------------------------
treatment_palette = sns.color_palette(PALETTE_NAME, n_colors=PALETTE_SIZE)
treatment_categories = None

# -----------------------------
# Helper Function: Custom HVG Plot
# -----------------------------
def plot_highly_variable_genes(adata, save_path="figures/04_hvg_colored.png"):
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
adata = sc.read(f"{RESULTS_DIR}/03_filtered_data.h5ad")
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

sc.pl.highly_variable_genes(adata, save="_04_hvg_plot.png", show=False)
plot_highly_variable_genes(adata, save_path=f"{FIGURE_DIR}/04_hvg_colored.png")
logging.info("Saved HVG plots")

# -----------------------------
# Scaling, PCA, and UMAP
# -----------------------------
sc.pp.scale(adata, max_value=SCALE_MAX)
logging.info(f"Scaled the data to unit variance with max value {SCALE_MAX}")

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color="treatment", palette=treatment_colors, save="_04_pca_treatment.png", show=False)
logging.info("Performed PCA and saved plot colored by treatment")

sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
sc.tl.umap(adata)
sc.pl.umap(adata, color="treatment", palette=treatment_colors, save="_04_umap_treatment.png", show=False)
logging.info(f"Computed neighbors (n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}) and UMAP, saved plot colored by treatment")

# -----------------------------
# Save Processed Data
# -----------------------------
adata.write(f"{RESULTS_DIR}/04_norm_dr.h5ad")
logging.info("Saved normalized and reduced AnnData object")

# -----------------------------
# Generate HTML Report
# -----------------------------
report_path = Path(f"{RESULTS_DIR}/04_report.html")
figures = sorted(Path(FIGURE_DIR).glob("*04_*.png"))

with open(report_path, "w") as f:
    f.write("<html><head><title>Step 4: Normalization & Dimensionality Reduction</title></head><body>\n")
    f.write("<h1>Step 4 Report: Normalization, HVG, PCA, UMAP</h1>\n")
    for fig in figures:
        f.write(f"<h3>{fig.name}</h3>\n")
        f.write(f'<img src="../{fig}" width="700"><br><br>\n')
    f.write("</body></html>")

logging.info(f"Generated HTML report at {report_path}")
print("✅ Step 4 complete: Normalization, HVG selection, PCA & UMAP done and saved.")
