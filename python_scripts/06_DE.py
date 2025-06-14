# -----------------------------
# Step 6: Differential Expression Analysis + Report
# -----------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pathlib import Path
import datetime
import os
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# User-adjustable parameters
# -----------------------------
RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("06_log.txt")
CSV_DIR = Path(f"{RESULTS_DIR}/DE_csv")
CSV_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = RESULTS_DIR / "05_clustered_data.h5ad"
OUTPUT_FILE = RESULTS_DIR / "06_de_data.h5ad"
HTML_REPORT = RESULTS_DIR / "06_report.html"

FIG_NAME1 = "06_DE"
FIG_NAME2 = "06_volcano"

CUTOFF_PVALS = 0.05
CUTOFF_DF = 1
RANK_GENES_METHOD = "wilcoxon"  # Marker gene ranking method
GROUPINGS = ['treatment', 'leiden']  # Grouping columns for DE
# -----------------------------


logging.info("Step 6 started: Differential Expression Analysis")

# -----------------------------
# Load data
# -----------------------------
adata = sc.read(INPUT_FILE)
logging.info(f"Loaded clustered data: {adata.n_obs} cells × {adata.n_vars} genes")

# -----------------------------
# Differential Expression
# -----------------------------
def plot_volcano(df, title, save_path):
    df['pvals_adj'] = df['pvals_adj'].fillna(1)
    df['log10_padj'] = -np.log10(df['pvals_adj'].replace(0, 1e-300))
    df['significant'] = (df['pvals_adj'] < CUTOFF_PVALS) & (df['logfoldchanges'].abs() > CUTOFF_DF)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x='logfoldchanges',
        y='log10_padj',
        hue='significant',
        palette={True: 'red', False: 'grey'},
        legend=False,
        alpha=0.7
    )
    plt.axhline(-np.log10(CUTOFF_PVALS), linestyle='--', color='blue')
    plt.axvline(1, linestyle='--', color='blue')
    plt.axvline(-1, linestyle='--', color='blue')
    plt.title(title)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-Value')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()

for groupby_col in GROUPINGS:
    logging.info(f"Running DE analysis grouped by '{groupby_col}'")
    key_added = f"rank_genes_{groupby_col}"
    sc.tl.rank_genes_groups(adata, groupby=groupby_col, method=RANK_GENES_METHOD, key_added=key_added, pts=True)

    unique_groups = sorted(adata.obs[groupby_col].unique())
    for grp in unique_groups:
        df = sc.get.rank_genes_groups_df(adata, group=grp, key=key_added)
        csv_path = CSV_DIR / f"{FIG_NAME1}_{groupby_col}_{grp}.csv"
        df.to_csv(csv_path, index=False)
        logging.info(f"Saved DE CSV for {groupby_col}={grp} at {csv_path}")

        fig_path = FIGURE_DIR / f"{FIG_NAME2}_{groupby_col}_{grp}.png"
        plot_volcano(df, f"Volcano Plot: {groupby_col}={grp}", fig_path)
        logging.info(f"Saved volcano plot for {groupby_col}={grp} at {fig_path}")

# -----------------------------
# Generate HTML Report
# -----------------------------
with open(HTML_REPORT, "w") as f:
    f.write("<html><head><title>Step 6: Differential Expression Analysis</title>\n")
    f.write("""
    <style>
        body { font-family: Arial, sans-serif; padding: 20px; }
        h1, h2, h3 { color: #2c3e50; }
        img { border: 1px solid #ccc; padding: 10px; margin-bottom: 20px; width: 600px; }
        a { text-decoration: none; color: #2980b9; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { padding: 8px 12px; border: 1px solid #ddd; }
    </style>
    </head><body>
    """)
    f.write("<h1>Step 6 Report: Differential Expression Analysis</h1>\n")
    f.write(f"<p><b>Date:</b> {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")

    for groupby_col in GROUPINGS:
        f.write(f"<h2>DE Results by {groupby_col.capitalize()}</h2>\n")
        for grp in sorted(adata.obs[groupby_col].unique()):
            csv_name = f"DE_csv/06_DE_{groupby_col}_{grp}.csv"
            fig_name = f"../figures/06_volcano_{groupby_col}_{grp}.png"
            df_preview = pd.read_csv(CSV_DIR / f"06_DE_{groupby_col}_{grp}.csv").head(10)
            f.write(f"<h3>{groupby_col.capitalize()}: {grp}</h3>\n")
            f.write(f"<p>Top 10 DE genes:</p>\n")
            f.write(df_preview.to_html(index=False, classes="table table-striped", border=0))
            f.write(f"<p>Full CSV: <a href=\"{csv_name}\">{csv_name}</a></p>\n")
            f.write(f"<img src=\"{fig_name}\" alt=\"Volcano plot {groupby_col} {grp}\"><br>\n")

    f.write("</body></html>\n")
logging.info(f"Generated HTML report at {HTML_REPORT}")

# -----------------------------
# Save updated AnnData
# -----------------------------
adata.write(OUTPUT_FILE)
logging.info(f"Saved DE results to {OUTPUT_FILE}")


print("✅ Step 6 complete: Differential expression analysis done and saved.")
