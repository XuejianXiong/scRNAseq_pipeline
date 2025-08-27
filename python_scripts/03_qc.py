# -----------------------------
# Step 3: Quality Control
# -----------------------------

import scanpy as sc
import pandas as pd
import numpy as np
from logzero import logger
from pathlib import Path
import sys
import datetime
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# User-Adjustable Parameters
# -----------------------------
#PROJ_NAME = "cropseq" 
PROJ_NAME = "retina"  

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("03_log.txt", PROJ_NAME)
INPUT_FILE = RESULTS_DIR / "02_merged_data.h5ad"
OUTPUT_FILE = RESULTS_DIR / "03_filtered_data.h5ad"
HTML_REPORT = RESULTS_DIR / "03_report.html"

MIN_GENES = 200
MAX_PCT_MT = 10
MIN_CELLS_PER_GENE = 3
    
# -----------------------------

# Control Scanpy output (independent of logzero)
sc.settings.verbosity = 1         # minimal output, no "saving figure..." warnings
sc.settings.logfile = LOG_FILE    # redirect scanpy logs to your log file

logger.info("Step 3 started: Quality Control")

# -----------------------------
# Load merged dataset from step 2
# -----------------------------
adata = sc.read(INPUT_FILE)
logger.info(f"Loaded merged data: {adata.n_obs} cells × {adata.n_vars} genes")

if 'treatment' in adata.obs.columns:
    GROUP_COLS = ['treatment']
elif 'sample' in adata.obs.columns:
    GROUP_COLS = ['sample']
else:
    GROUP_COLS = []
    logger.warning("No 'treatment' or 'sample' column found in adata.obs. Exiting.")
    sys.exit(1)


# -----------------------------
# Annotate mitochondrial genes and calculate QC metrics
# -----------------------------
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# -----------------------------
# Save QC plots before filtering
# -----------------------------
for group in GROUP_COLS:
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        groupby=group,
        jitter=0.4,
        rotation=45,
        save=f"_qc_violin_by_{group}.png",
        show=False
    )
    logger.info(f"Saved violin plot grouped by {group} before filtering")

    sc.pl.scatter(
        adata,
        x='total_counts',
        y='pct_counts_mt',
        color=group,
        save=f"_qc_scatter_mt_by_{group}_before.png",
        show=False
    )
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color=group,
        save=f"_qc_scatter_genes_by_{group}_before.png",
        show=False
    )
    logger.info(f"Saved scatter plots colored by {group} before filtering")

pre_shape = adata.shape

# -----------------------------
# Apply filters
# -----------------------------
sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
adata = adata[adata.obs['pct_counts_mt'] < MAX_PCT_MT, :]

post_shape = adata.shape
logger.info(f"Filtering complete — Before: {pre_shape}, After: {post_shape}")
logger.info(f"Applied thresholds — min_genes: {MIN_GENES}, min_cells_per_gene: {MIN_CELLS_PER_GENE}, max_pct_mito: {MAX_PCT_MT}%")

# -----------------------------
# Save QC plots after filtering
# -----------------------------
for group in GROUP_COLS:
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        groupby=group,
        jitter=0.4,
        rotation=45,
        save=f"_qc_violin_by_{group}_after.png",
        show=False
    )
    logger.info(f"Saved violin plot grouped by {group} after filtering")

    sc.pl.scatter(
        adata,
        x='total_counts',
        y='pct_counts_mt',
        color=group,
        save=f"_qc_scatter_mt_by_{group}_after.png",
        show=False
    )
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color=group,
        save=f"_qc_scatter_genes_by_{group}_after.png",
        show=False
    )
    logger.info(f"Saved scatter plots colored by {group} after filtering")

# -----------------------------
# Save filtered data
# -----------------------------
adata.write(OUTPUT_FILE)
logger.info(f"Filtered data saved to {OUTPUT_FILE}")

# -----------------------------
# Generate HTML report with embedded images
# -----------------------------
# Use absolute paths in the report HTML
FIGURE_DIR_ABS = FIGURE_DIR.resolve()

with open(HTML_REPORT, "w") as f:
    f.write(f"""<html><head><title>QC Report</title></head><body>
<h1>Single-Cell QC Report</h1>
<p><b>Date:</b> {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
<p><b>Initial shape:</b> {pre_shape[0]} cells × {pre_shape[1]} genes</p>
<p><b>Post-filtering shape:</b> {post_shape[0]} cells × {post_shape[1]} genes</p>
<p><b>QC thresholds:</b> min_genes = {MIN_GENES}, min_cells_per_gene = {MIN_CELLS_PER_GENE}, max_pct_mito = {MAX_PCT_MT}%</p>

<h2>QC Metrics Before Filtering</h2>
<h3>Violin plots grouped by metadata</h3>
""")

    for group in GROUP_COLS:
        f.write(f'<img src="{FIGURE_DIR_ABS}/violin_qc_violin_by_{group}.png" width="600"/><br>\n')

    f.write("<h3>Scatter plots colored by metadata</h3>\n")

    for group in GROUP_COLS:
        f.write(f'<img src="{FIGURE_DIR_ABS}/scatter_qc_scatter_mt_by_{group}_before.png" width="600"/><br>\n')
        f.write(f'<img src="{FIGURE_DIR_ABS}/scatter_qc_scatter_genes_by_{group}_before.png" width="600"/><br>\n')

    f.write("""
<h2>QC Metrics After Filtering</h2>
<h3>Violin plots grouped by metadata</h3>
""")

    for group in GROUP_COLS:
        f.write(f'<img src="{FIGURE_DIR_ABS}/violin_qc_violin_by_{group}_after.png" width="600"/><br>\n')

    f.write("<h3>Scatter plots colored by metadata</h3>\n")

    for group in GROUP_COLS:
        f.write(f'<img src="{FIGURE_DIR_ABS}/scatter_qc_scatter_mt_by_{group}_after.png" width="600"/><br>\n')
        f.write(f'<img src="{FIGURE_DIR_ABS}/scatter_qc_scatter_genes_by_{group}_after.png" width="600"/><br>\n')

    f.write("</body></html>\n")
logger.info(f"HTML QC report generated at: {HTML_REPORT}")


print("✅ Step 3 complete: Quality control applied, filtered data saved.")
