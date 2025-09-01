#!/usr/bin/env python3
"""
Step 5: Clustering and Marker Gene Analysis for scRNA-seq data.

This step:
- Performs clustering using the Leiden algorithm.
- Computes neighborhood graph for cells.
- Identifies marker genes per cluster.
- Generates plots (UMAP colored by cluster, marker gene heatmaps).
- Creates an HTML summary report.
"""

import scanpy as sc
import pandas as pd
from logzero import logger
from pathlib import Path
import yaml
import sys
import os
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action="ignore", category=FutureWarning)


# -----------------------------
# Load Config and Parameters
# -----------------------------
if len(sys.argv) < 3:
    print("Usage: step5_clustering.py <project_name> <config.yaml>")
    sys.exit(1)

PROJ_NAME = sys.argv[1]
CONFIG_FILE = Path(sys.argv[2])

if not CONFIG_FILE.exists():
    print(f"Config file not found: {CONFIG_FILE}")
    sys.exit(1)

with open(CONFIG_FILE, "r") as f:
    cfg = yaml.safe_load(f)

# Enforce strict clustering parameters
try:
    CLUSTER_PARAMS = cfg["clustering"]
    LEIDEN_RESOLUTION = CLUSTER_PARAMS["leiden_resolution"]
    N_NEIGHBORS = CLUSTER_PARAMS["n_neighbors"]
    N_PCS = CLUSTER_PARAMS["n_pcs"]
    RANK_GENES_METHOD = CLUSTER_PARAMS["rank_genes_method"]
    TOP_N_GENES = CLUSTER_PARAMS["top_n_genes"]
except KeyError as e:
    print(f"Missing clustering parameter in config.yaml: {e}")
    sys.exit(1)

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("05_log.txt", PROJ_NAME)
MARKER_DIR = RESULTS_DIR / "markers"
MARKER_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = RESULTS_DIR / "04_norm_dr_data.h5ad"
CSV_FILE = RESULTS_DIR / "05_top_marker_genes_summary.csv"
OUTPUT_FILE = RESULTS_DIR / "05_clustered_data.h5ad"
HTML_REPORT = RESULTS_DIR / "05_report.html"

sc.settings.verbosity = 1
sc.settings.logfile = LOG_FILE
sc.settings.figdir = FIGURE_DIR

logger.info("Step 5 started: Clustering and Marker Gene Analysis")


# -----------------------------
# Helper Functions
# -----------------------------
def identify_and_save_marker_genes(adata, marker_dir, csv_file, rank_method, top_n):
    """
    Identify marker genes for each cluster, save per-cluster CSV files,
    and save a summary CSV of top N marker genes per cluster.
    """
    logger.info(f"Computing marker genes per cluster using {rank_method} method")
    sc.tl.rank_genes_groups(adata, groupby="cluster", method=rank_method)

    clusters = (
        adata.obs["cluster"].cat.categories
        if hasattr(adata.obs["cluster"], "cat")
        else sorted(adata.obs["cluster"].unique())
    )

    summary_rows = []
    for cluster in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df["cluster"] = cluster
        df = df[["cluster", "names", "logfoldchanges", "pvals_adj"]]
        df.columns = ["cluster", "gene", "logfoldchanges", "pvals_adj"]
        out_file = marker_dir / f"marker_genes_cluster_{cluster}.csv"
        df.to_csv(out_file, index=False)
        logger.info(f"Saved marker genes for cluster {cluster} with {len(df)} genes")
        summary_rows.append(df.head(top_n))

    summary_df = pd.concat(summary_rows, ignore_index=True)
    summary_df.to_csv(csv_file, index=False)
    logger.info(f"Saved top {top_n} marker genes summary for all clusters.")
    return clusters


def generate_html_report(html_file, figures, clusters, top_n):
    with open(html_file, "w") as f:
        f.write(
            "<html><head><title>Step 5: Clustering & Marker Gene Analysis</title>\n"
        )
        f.write(
            """
        <style>
            body { font-family: Arial, sans-serif; padding: 20px; }
            h1, h2 { color: #2c3e50; }
            img { border: 1px solid #ccc; padding: 10px; margin-bottom: 20px; width: 700px; }
            a { text-decoration: none; color: #2980b9; }
        </style>
        </head><body>
        """
        )
        f.write("<h1>Step 5 Report: Leiden Clustering and Marker Genes</h1>\n")

        f.write("<h2>UMAP Plots</h2>\n")
        for fig in figures:
            rel_path = os.path.relpath(fig, start=html_file.parent)
            f.write(f'<img src="{rel_path}" alt="{fig.name}"><br>\n')

        f.write("<h2>Marker Genes per Cluster</h2>\n<ul>\n")
        for cluster in clusters:
            f.write(
                f'<li><a href="markers/marker_genes_cluster_{cluster}.csv">Marker genes for cluster {cluster}</a></li>\n'
            )
        f.write("</ul>\n")

        f.write(
            f'<h2>Summary</h2>\n<p><a href="05_top_marker_genes_summary.csv">Top {top_n} marker genes summary (CSV)</a></p>\n'
        )
        f.write("</body></html>")
    logger.info(f"Generated HTML report at {html_file}")


# -----------------------------
# Main Pipeline
# -----------------------------
def main() -> None:
    adata = sc.read(INPUT_FILE)
    logger.info(f"Loaded normalized data: {adata.n_obs} cells, {adata.n_vars} genes")

    # Recompute neighbors
    logger.info(f"Computing neighbors with n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)

    # Leiden clustering
    logger.info(f"Running Leiden clustering with resolution={LEIDEN_RESOLUTION}")
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION)
    adata.obs["cluster"] = adata.obs["leiden"]

    cluster_sizes = adata.obs["cluster"].value_counts().sort_index()
    logger.info(f"Cluster sizes:\n{cluster_sizes.to_string()}")

    # UMAP plots
    logger.info("Generating UMAP plots colored by clusters and treatment/sample")
    sc.pl.umap(adata, color="cluster", save="_05_umap_clusters.png", show=False)
    if "treatment" in adata.obs.columns:
        sc.pl.umap(adata, color="treatment", save="_05_umap_treatment.png", show=False)
    elif "sample" in adata.obs.columns:
        sc.pl.umap(adata, color="sample", save="_05_umap_sample.png", show=False)
    else:
        logger.warning("No 'treatment' or 'sample' metadata found, skipping plot.")

    # Marker gene identification + summary
    clusters = identify_and_save_marker_genes(
        adata, MARKER_DIR, CSV_FILE, RANK_GENES_METHOD, TOP_N_GENES
    )

    # Generate HTML report
    figures = sorted(Path(sc.settings.figdir).glob("*05_*.png"))
    generate_html_report(HTML_REPORT, figures, clusters, TOP_N_GENES)

    # Save clustered AnnData object
    adata.write(OUTPUT_FILE)
    logger.info(f"Saved clustered AnnData object to {OUTPUT_FILE}")

    logger.info(
        "✅ Step 5 complete: Clustering and marker gene analysis done and saved."
    )


if __name__ == "__main__":
    main()
