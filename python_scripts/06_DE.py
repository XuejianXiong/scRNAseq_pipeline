#!/usr/bin/env python3
"""
Step 6: Differential Expression Analysis for scRNA-seq data.

This step:
- Runs differential expression (DE) analysis on clusters and sample/treatment groups.
- Saves DE results as CSV files and generates volcano plots.
- Creates an HTML summary report of top DE genes.
"""

import sys
import warnings
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
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
    print("Usage: step6_de.py <project_name> <config.yaml>")
    sys.exit(1)

PROJ_NAME = sys.argv[1]
CONFIG_FILE = Path(sys.argv[2])

if not CONFIG_FILE.exists():
    print(f"Config file not found: {CONFIG_FILE}")
    sys.exit(1)

with open(CONFIG_FILE, "r") as f:
    cfg = yaml.safe_load(f)

# Enforce strict DE thresholds
try:
    DE_PARAMS = cfg["de"]
    CUTOFF_PVALS = DE_PARAMS["cutoff_pvals"]
    CUTOFF_LFC = DE_PARAMS["cutoff_lfc"]
    RANK_GENES_METHOD = DE_PARAMS["rank_genes_method"]
except KeyError as e:
    print(f"Missing DE parameter in config.yaml: {e}")
    sys.exit(1)

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("06_log.txt", PROJ_NAME)
INPUT_FILE = RESULTS_DIR / "05_clustered_data.h5ad"
OUTPUT_FILE = RESULTS_DIR / "06_de_data.h5ad"
HTML_REPORT = RESULTS_DIR / "06_report.html"
CSV_DIR = RESULTS_DIR / "DE_csv"
CSV_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1  # minimal console output
sc.settings.logfile = LOG_FILE  # direct Scanpy logs to log file
sc.settings.figdir = FIGURE_DIR  # save plots automatically here

logger.info("Step 6 started: Differential Expression Analysis")


# -----------------------------
# Helper Functions
# -----------------------------
def plot_volcano(df: pd.DataFrame, title: str, save_path: Path) -> None:
    """Generate a volcano plot from DE results."""
    df = df.copy()
    df["pvals_adj"] = df["pvals_adj"].fillna(1)
    df["log10_padj"] = -np.log10(df["pvals_adj"].replace(0, 1e-300))
    df["significant"] = (df["pvals_adj"] < CUTOFF_PVALS) & (
        df["logfoldchanges"].abs() > CUTOFF_LFC
    )

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="logfoldchanges",
        y="log10_padj",
        hue="significant",
        palette={True: "red", False: "grey"},
        legend=False,
        alpha=0.7,
    )
    plt.axhline(-np.log10(CUTOFF_PVALS), linestyle="--", color="blue")
    plt.axvline(CUTOFF_LFC, linestyle="--", color="blue")
    plt.axvline(-CUTOFF_LFC, linestyle="--", color="blue")
    plt.title(title)
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted P-Value")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


def run_de_analysis(adata: sc.AnnData, groupby_col: str) -> None:
    """Perform DE analysis per grouping and save CSV + volcano plots."""
    key_added = f"rank_genes_{groupby_col}"
    logger.info(
        f"Running DE analysis grouped by '{groupby_col}' using {RANK_GENES_METHOD}"
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby_col,
        method=RANK_GENES_METHOD,
        key_added=key_added,
        pts=True,
    )

    for grp in sorted(adata.obs[groupby_col].unique()):
        df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)
        csv_path = CSV_DIR / f"06_DE_{groupby_col}_{grp}.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"Saved DE CSV for {groupby_col}={grp} at {csv_path}")

        fig_path = FIGURE_DIR / f"06_volcano_{groupby_col}_{grp}.png"
        plot_volcano(df, f"Volcano Plot: {groupby_col}={grp}", fig_path)
        logger.info(f"Saved volcano plot for {groupby_col}={grp} at {fig_path}")


def generate_html_report(fig_dir: Path, output_file: Path, groupings: list) -> None:
    """Generate an HTML report embedding volcano plots and DE previews."""
    with open(output_file, "w") as f:
        f.write("<html><head><title>Step 6: Differential Expression</title>\n")
        f.write(
            "<style>body{font-family:Arial;padding:20px;}img{width:600px;margin:10px;}</style></head><body>\n"
        )
        f.write("<h1>Step 6 Report: Differential Expression Analysis</h1>\n")
        f.write(
            f"<p><b>Date:</b> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n"
        )

        fig_dir_abs = fig_dir.resolve()
        for col in groupings:
            f.write(f"<h2>DE Results by {col.capitalize()}</h2>\n")
            for grp in sorted(adata.obs[col].unique()):
                csv_name = f"DE_csv/06_DE_{col}_{grp}.csv"
                fig_name = fig_dir_abs / f"06_volcano_{col}_{grp}.png"
                df_preview = pd.read_csv(CSV_DIR / f"06_DE_{col}_{grp}.csv").head(10)

                f.write(f"<h3>{col.capitalize()}: {grp}</h3>\n")
                f.write("<p>Top 10 DE genes:</p>\n")
                f.write(df_preview.to_html(index=False, border=0))
                f.write(f"<p>Full CSV: <a href='{csv_name}'>{csv_name}</a></p>\n")
                f.write(f"<img src='{fig_name}' alt='Volcano plot {col} {grp}'><br>\n")
        f.write("</body></html>\n")
    logger.info(f"Generated HTML report at {output_file}")


# -----------------------------
# Main Pipeline
# -----------------------------
def main() -> None:
    global adata  # Used inside generate_html_report
    adata = sc.read(INPUT_FILE)
    logger.info(f"Loaded clustered data: {adata.n_obs} cells × {adata.n_vars} genes")

    # Determine grouping columns
    if "treatment" in adata.obs.columns:
        groupings = ["treatment", "leiden"]
    elif "sample" in adata.obs.columns:
        groupings = ["sample", "leiden"]
    else:
        logger.error("No 'treatment' or 'sample' columns found — cannot run DE.")
        sys.exit(1)

    # Run DE analysis for each grouping
    for col in groupings:
        run_de_analysis(adata, col)

    # Save updated AnnData object
    adata.write(OUTPUT_FILE)
    logger.info(f"Saved DE results to {OUTPUT_FILE}")

    # Generate HTML report
    generate_html_report(FIGURE_DIR, HTML_REPORT, groupings)

    logger.info("✅ Step 6 complete: Differential expression analysis done and saved.")


if __name__ == "__main__":
    main()
