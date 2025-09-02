#!/usr/bin/env python3
"""
Step 4: Normalization and Dimensionality Reduction for scRNA-seq data.

This step:
- Normalizes counts per cell.
- Applies log-transformation.
- Selects highly variable genes (HVGs).
- Performs scaling, PCA, and UMAP.
- Generates plots and an HTML summary report.
"""

import sys
import warnings
from pathlib import Path
from typing import Optional

import yaml
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from logzero import logger

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action="ignore", category=FutureWarning)

# -----------------------------
# Load Config and Parameters
# -----------------------------
if len(sys.argv) < 3:
    print("Usage: step4_norm_dr.py <project_name> <config.yaml>")
    sys.exit(1)

PROJ_NAME = sys.argv[1]
CONFIG_FILE = Path(sys.argv[2])

if not CONFIG_FILE.exists():
    print(f"Config file not found: {CONFIG_FILE}")
    sys.exit(1)

with open(CONFIG_FILE, "r") as f:
    cfg = yaml.safe_load(f)

# Enforce strict normalization parameters
try:
    NORM_PARAMS = cfg["normalization"]
    N_TOP_HVG = NORM_PARAMS["n_top_hvg"]
    SCALE_MAX = NORM_PARAMS["scale_max"]
    N_PCS = NORM_PARAMS["n_pcs"]
    N_NEIGHBORS = NORM_PARAMS["n_neighbors"]
except KeyError as e:
    print(f"Missing normalization parameter in config.yaml: {e}")
    sys.exit(1)

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("04_log.txt", PROJ_NAME)
INPUT_FILE = RESULTS_DIR / "03_filtered_data.h5ad"
OUTPUT_FILE = RESULTS_DIR / "04_norm_dr_data.h5ad"
HTML_REPORT = RESULTS_DIR / "04_report.html"

PALETTE_NAME = "Set2"
PALETTE_SIZE = 10
PALETTE = sns.color_palette(PALETTE_NAME, n_colors=PALETTE_SIZE)

sc.settings.verbosity = 1  # minimal console output
sc.settings.logfile = LOG_FILE  # direct Scanpy logs to log file
sc.settings.figdir = FIGURE_DIR  # save plots automatically here

logger.info("Step 4 started: Normalization and Dimensionality Reduction")


# -----------------------------
# Helper Functions
# -----------------------------
def plot_highly_variable_genes(adata: sc.AnnData, save_path: Path) -> None:
    """Custom HVG plot (mean vs. dispersion)."""
    if "highly_variable" not in adata.var:
        raise ValueError("HVGs not computed before plotting.")

    means = adata.var["means"]
    dispersions = adata.var["dispersions"]
    hvg_mask = adata.var["highly_variable"]

    plt.figure(figsize=(6, 4))
    plt.scatter(
        means[~hvg_mask],
        dispersions[~hvg_mask],
        c="lightgray",
        s=10,
        label="Not HVG",
        alpha=0.5,
    )
    plt.scatter(
        means[hvg_mask], dispersions[hvg_mask], c="red", s=10, label="HVG", alpha=0.8
    )
    plt.xlabel("Mean expression")
    plt.ylabel("Dispersion")
    plt.title("Highly Variable Genes")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


def detect_grouping_column(adata: sc.AnnData) -> Optional[str]:
    """Determine the best obs column for group coloring."""
    for col in ["treatment", "sample"]:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype("category")
            return col
    return None


def assign_group_colors(adata: sc.AnnData, group_col: str) -> dict:
    """Assign a color palette to group categories."""
    categories = adata.obs[group_col].cat.categories
    colors = dict(zip(categories, PALETTE[: len(categories)]))
    adata.uns[f"{group_col}_colors"] = [colors[c] for c in categories]
    return colors


def generate_html_report(fig_dir: Path, output_file: Path) -> None:
    """Generate an HTML report embedding Step 4 plots."""
    figs = sorted(fig_dir.resolve().glob("*04_*.png"))
    with open(output_file, "w") as f:
        f.write("<html><head><title>Step 4: Normalization & DR</title></head><body>\n")
        f.write("<h1>Step 4 Report: Normalization, HVG, PCA, UMAP</h1>\n")
        for fig in figs:
            f.write(f"<h3>{fig.name}</h3>\n")
            f.write(f'<img src="{str(fig)}" width="700"><br><br>\n')
        f.write("</body></html>")


# -----------------------------
# Main Pipeline
# -----------------------------
def main() -> None:
    # Load data
    adata = sc.read(INPUT_FILE)
    logger.info(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")

    # Determine grouping column for plots
    group_col = detect_grouping_column(adata)
    group_colors = assign_group_colors(adata, group_col) if group_col else None
    if group_col:
        logger.info(
            f"Using '{group_col}' as grouping column with colors: {group_colors}"
        )
    else:
        logger.warning(
            "No 'treatment' or 'sample' column found — plots will not be colored by group."
        )

    # Normalization & HVG selection
    sc.pp.normalize_total(adata, target_sum=1e4)
    logger.info("Normalized total counts per cell to 1e4")

    sc.pp.log1p(adata)
    adata.raw = adata
    logger.info("Applied log1p transformation")

    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_TOP_HVG, flavor="seurat", subset=True
    )
    logger.info(
        f"Selected {np.sum(adata.var['highly_variable'])} highly variable genes"
    )

    plot_highly_variable_genes(adata, FIGURE_DIR / "04_hvg_custom.png")
    sc.pl.highly_variable_genes(adata, save="_04_hvg_scanpy.png", show=False)
    logger.info("Saved HVG plots")

    # Scaling, PCA, and UMAP
    sc.pp.scale(adata, max_value=SCALE_MAX)
    logger.info(f"Scaled data to unit variance with max value {SCALE_MAX}")

    sc.tl.pca(adata, svd_solver="arpack")
    sc.pl.pca(
        adata,
        color=group_col if group_col else None,
        palette=group_colors if group_colors else None,
        save="_04_pca.png",
        show=False,
    )
    logger.info("Performed PCA and saved plot")

    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color=group_col if group_col else None,
        palette=group_colors if group_colors else None,
        save="_04_umap.png",
        show=False,
    )
    logger.info(
        f"Computed neighbors (n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}) and UMAP, saved plot"
    )

    # Save processed data
    adata.write(OUTPUT_FILE)
    logger.info(f"Saved AnnData to {OUTPUT_FILE}")

    # Generate HTML report
    generate_html_report(FIGURE_DIR, HTML_REPORT)
    logger.info(f"Generated HTML report at {HTML_REPORT}")

    logger.info(
        "✅ Step 4 complete: Normalization, HVG selection, PCA & UMAP done and saved."
    )


if __name__ == "__main__":
    main()
