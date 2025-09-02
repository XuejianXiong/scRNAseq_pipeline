#!/usr/bin/env python3
"""
Step 3: Quality Control for single-cell RNA-seq data.

- Annotate QC metrics (gene counts, mitochondrial percentage).
- Generate violin and scatter plots grouped by metadata.
- Filter cells and genes based on quality thresholds.
- Save filtered data and generate an HTML QC report.
"""

import sys
import yaml
import datetime
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from logzero import logger
from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action="ignore", category=FutureWarning)


# -----------------------------
# Load Config and Parameters
# -----------------------------
if len(sys.argv) < 3:
    print("Usage: step3_qc.py <project_name> <config.yaml>")
    sys.exit(1)

PROJ_NAME = sys.argv[1]
CONFIG_FILE = Path(sys.argv[2])

if not CONFIG_FILE.exists():
    print(f"Config file not found: {CONFIG_FILE}")
    sys.exit(1)

with open(CONFIG_FILE, "r") as f:
    cfg = yaml.safe_load(f)

# Enforce strict QC parameters
try:
    QC_PARAMS = cfg["qc"]
    MIN_GENES = QC_PARAMS["min_genes"]
    MAX_PCT_MT = QC_PARAMS["max_pct_mt"]
    MIN_CELLS_PER_GENE = QC_PARAMS["min_cells_per_gene"]
except KeyError as e:
    print(f"Missing QC parameter in config.yaml: {e}")
    sys.exit(1)

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("03_log.txt", PROJ_NAME)
INPUT_FILE = RESULTS_DIR / "02_merged_data.h5ad"
OUTPUT_FILE = RESULTS_DIR / "03_filtered_data.h5ad"
HTML_REPORT = RESULTS_DIR / "03_report.html"

sc.settings.verbosity = 1  # minimal console output
sc.settings.logfile = LOG_FILE  # direct Scanpy logs to log file
sc.settings.figdir = FIGURE_DIR  # save plots automatically here

logger.info("Step 3 started: Quality Control")


# -----------------------------
# Helper Functions
# -----------------------------
def load_data(input_file: Path):
    """Load AnnData object from previous step."""
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)
    adata = sc.read(input_file)
    logger.info(f"Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def determine_grouping(adata):
    """Choose grouping column for QC plots."""
    if "treatment" in adata.obs.columns:
        return ["treatment"]
    elif "sample" in adata.obs.columns:
        return ["sample"]
    logger.error("No 'treatment' or 'sample' column found in adata.obs.")
    sys.exit(1)


def annotate_qc(adata):
    """Add mitochondrial gene flag and compute QC metrics."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


def save_qc_plots(adata, group_cols, tag):
    """Save QC violin and scatter plots for each grouping column."""
    for group in group_cols:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            groupby=group,
            jitter=0.4,
            rotation=45,
            save=f"_qc_violin_by_{group}_{tag}.png",
            show=False,
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="pct_counts_mt",
            color=group,
            save=f"_qc_scatter_mt_by_{group}_{tag}.png",
            show=False,
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="n_genes_by_counts",
            color=group,
            save=f"_qc_scatter_genes_by_{group}_{tag}.png",
            show=False,
        )
        logger.info(f"Saved QC plots grouped by {group} ({tag} filtering)")


def filter_data(adata):
    """Apply cell/gene QC filters and return filtered AnnData."""
    pre_shape = adata.shape
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MT, :]
    post_shape = adata.shape

    logger.info(f"Filtering complete — Before: {pre_shape}, After: {post_shape}")
    logger.info(
        f"Applied thresholds — min_genes={MIN_GENES}, "
        f"min_cells_per_gene={MIN_CELLS_PER_GENE}, max_pct_mito={MAX_PCT_MT}%"
    )
    return adata, pre_shape, post_shape


def save_html_report(html_file, figure_dir, pre_shape, post_shape, group_cols):
    """Generate a static HTML report embedding QC plots."""
    figure_dir = figure_dir.resolve()
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(html_file, "w") as f:
        f.write(f"<html><head><title>QC Report</title></head><body>\n")
        f.write("<h1>Single-Cell QC Report</h1>\n")
        f.write(f"<p><b>Date:</b> {now}</p>\n")
        f.write(
            f"<p><b>Initial shape:</b> {pre_shape[0]} cells × {pre_shape[1]} genes</p>\n"
        )
        f.write(
            f"<p><b>Post-filtering shape:</b> {post_shape[0]} cells × {post_shape[1]} genes</p>\n"
        )
        f.write(
            f"<p><b>QC thresholds:</b> min_genes={MIN_GENES}, "
            f"min_cells_per_gene={MIN_CELLS_PER_GENE}, "
            f"max_pct_mito={MAX_PCT_MT}%</p>\n"
        )

        f.write("<h2>QC Metrics Before Filtering</h2>\n")
        for group in group_cols:
            f.write(
                f'<img src="{figure_dir}/violin_qc_violin_by_{group}_before.png" width="600"><br>\n'
            )
            f.write(
                f'<img src="{figure_dir}/scatter_qc_scatter_mt_by_{group}_before.png" width="600"><br>\n'
            )
            f.write(
                f'<img src="{figure_dir}/scatter_qc_scatter_genes_by_{group}_before.png" width="600"><br>\n'
            )

        f.write("<h2>QC Metrics After Filtering</h2>\n")
        for group in group_cols:
            f.write(
                f'<img src="{figure_dir}/violin_qc_violin_by_{group}_after.png" width="600"><br>\n'
            )
            f.write(
                f'<img src="{figure_dir}/scatter_qc_scatter_mt_by_{group}_after.png" width="600"><br>\n'
            )
            f.write(
                f'<img src="{figure_dir}/scatter_qc_scatter_genes_by_{group}_after.png" width="600"><br>\n'
            )

        f.write("</body></html>\n")
    logger.info(f"HTML QC report generated at: {html_file}")


# -----------------------------
# Main Pipeline
# -----------------------------
def main():
    adata = load_data(INPUT_FILE)
    group_cols = determine_grouping(adata)

    annotate_qc(adata)
    save_qc_plots(adata, group_cols, tag="before")

    adata_filtered, pre_shape, post_shape = filter_data(adata)
    save_qc_plots(adata_filtered, group_cols, tag="after")

    adata_filtered.write(OUTPUT_FILE)
    logger.info(f"Filtered data saved to {OUTPUT_FILE}")

    save_html_report(HTML_REPORT, FIGURE_DIR, pre_shape, post_shape, group_cols)
    logger.info("Step 3 complete: QC applied, filtered data saved.")


if __name__ == "__main__":
    main()
