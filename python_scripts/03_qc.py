import scanpy as sc
import pandas as pd
import numpy as np
import os
import datetime

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# Step 3: Quality Control
# -----------------------------

# Define paths
input_file = "results/02_merged_data.h5ad"
output_file = "results/03_filtered_data.h5ad"
qc_plot_dir = "figures"
html_report = "results/03_report.html"
log_file = "results/03_log.txt"

# QC thresholds
min_genes = 200
max_pct_mt = 10
min_cells_per_gene = 3

# Create output directories
os.makedirs(qc_plot_dir, exist_ok=True)
sc.settings.figdir = qc_plot_dir
sc.settings.autoshow = False

# Initialize log
def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{timestamp}] {msg}"
    print(full_msg)
    with open(log_file, "a") as f:
        f.write(full_msg + "\n")

# Load merged dataset
adata = sc.read(input_file)
log(f"Loaded merged data: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Annotate mitochondrial genes and calculate QC metrics
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Metadata columns to use for grouping/coloring
group_cols = ['treatment']

# Save QC plots before filtering

# Violin plots grouped by metadata
for group in group_cols:
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        groupby=group,
        jitter=0.4,
        rotation=45,
        save=f"_qc_violin_by_{group}.png"
    )
    log(f"Saved violin plot grouped by {group} before filtering")

# Scatter plots colored by metadata
for color in group_cols:
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='pct_counts_mt',
        color=color,
        save=f"_qc_scatter_mt_by_{color}_before.png"
    )
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color=color,
        save=f"_qc_scatter_genes_by_{color}_before.png"
    )
    log(f"Saved scatter plots colored by {color} before filtering")

# Record shape before filtering
pre_shape = adata.shape

# Apply filters
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
adata = adata[adata.obs['pct_counts_mt'] < max_pct_mt, :]

# Record shape after filtering
post_shape = adata.shape
log(f"Filtering complete — Before: {pre_shape}, After: {post_shape}")
log(f"Applied thresholds — min_genes: {min_genes}, min_cells_per_gene: {min_cells_per_gene}, max_pct_mito: {max_pct_mt}%")

# Save QC plots after filtering

# Violin plots grouped by metadata after filtering
for group in group_cols:
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        groupby=group,
        jitter=0.4,
        rotation=45,
        save=f"_qc_violin_by_{group}_after.png"
    )
    log(f"Saved violin plot grouped by {group} after filtering")

# Scatter plots colored by metadata after filtering
for color in group_cols:
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='pct_counts_mt',
        color=color,
        save=f"_qc_scatter_mt_by_{color}_after.png"
    )
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color=color,
        save=f"_qc_scatter_genes_by_{color}_after.png"
    )
    log(f"Saved scatter plots colored by {color} after filtering")

# Save filtered data
adata.write(output_file)
log(f"Filtered data saved to {output_file}")

# Generate HTML report including grouped plots
with open(html_report, "w") as f:
    f.write(f"""
    <html><head><title>QC Report</title></head><body>
    <h1>Single-Cell QC Report</h1>
    <p><b>Date:</b> {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    <p><b>Initial shape:</b> {pre_shape[0]} cells × {pre_shape[1]} genes</p>
    <p><b>Post-filtering shape:</b> {post_shape[0]} cells × {post_shape[1]} genes</p>
    <p><b>QC thresholds:</b> min_genes = {min_genes}, min_cells_per_gene = {min_cells_per_gene}, max_pct_mito = {max_pct_mt}%</p>

    <h2>QC Metrics Before Filtering</h2>
    <h3>Violin plots grouped by metadata</h3>
    """)

    for group in group_cols:
        f.write(f'<img src="../{qc_plot_dir}/violin_qc_violin_by_{group}.png" width="600"/><br>\n')

    f.write("<h3>Scatter plots colored by metadata</h3>\n")

    for color in group_cols:
        f.write(f'<img src="../{qc_plot_dir}/scatter_qc_scatter_mt_by_{color}_before.png" width="600"/><br>\n')
        f.write(f'<img src="../{qc_plot_dir}/scatter_qc_scatter_genes_by_{color}_before.png" width="600"/><br>\n')

    f.write("""
    <h2>QC Metrics After Filtering</h2>
    <h3>Violin plots grouped by metadata</h3>
    """)

    for group in group_cols:
        f.write(f'<img src="../{qc_plot_dir}/violin_qc_violin_by_{group}_after.png" width="600"/><br>\n')

    f.write("<h3>Scatter plots colored by metadata</h3>\n")

    for color in group_cols:
        f.write(f'<img src="../{qc_plot_dir}/scatter_qc_scatter_mt_by_{color}_after.png" width="600"/><br>\n')
        f.write(f'<img src="../{qc_plot_dir}/scatter_qc_scatter_genes_by_{color}_after.png" width="600"/><br>\n')

    f.write("</body></html>\n")

log(f"HTML QC report generated at: {html_report}")

print("✅ Step 3 complete: Quality control applied, filtered data saved.")
