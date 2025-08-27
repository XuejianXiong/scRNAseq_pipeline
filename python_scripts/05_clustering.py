# -----------------------------
# Step 5: Clustering and Marker Gene Analysis
# -----------------------------

import scanpy as sc
import pandas as pd
from logzero import logger
from pathlib import Path
import os
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# User-adjustable parameters
# -----------------------------
#PROJ_NAME = "cropseq"  
PROJ_NAME = "retina"  

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("05_log.txt", PROJ_NAME)
MARKER_DIR = Path(f"{RESULTS_DIR}/markers")
MARKER_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = RESULTS_DIR / "04_norm_dr_data.h5ad"
CSV_FILE = RESULTS_DIR / "05_top_marker_genes_summary.csv"
OUTPUT_FILE = RESULTS_DIR / "05_clustered_data.h5ad"
HTML_REPORT = RESULTS_DIR / "05_report.html"

FIG_NAME1 = "05_umap_clusters.png"
FIG_NAME2 = "05_umap_treatment.png"

LEIDEN_RESOLUTION = 0.5        # Adjust clustering granularity (higher = more clusters)
N_NEIGHBORS = 15               # Number of neighbors to use for graph construction (adjust based on data)
N_PCS = 40                     # Number of PCs to use
RANK_GENES_METHOD = "wilcoxon"  # Marker gene ranking method
TOP_N_GENES = 10               # Number of top marker genes per cluster to save in summary
# -----------------------------


logger.info("Step 5 started: Clustering and Marker Gene Analysis")

# -----------------------------
# Load normalized and reduced data from Step 4
# -----------------------------
adata = sc.read(INPUT_FILE)
logger.info(f"Loaded normalized data: {adata.n_obs} cells, {adata.n_vars} genes")

# -----------------------------
# Recompute neighbors
# -----------------------------
logger.info(f"Computing neighbors with n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}")
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)

# -----------------------------
# Clustering with Leiden algorithm
# -----------------------------
logger.info(f"Running Leiden clustering with resolution={LEIDEN_RESOLUTION}")
sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION)

# Save clustering labels in obs
adata.obs["cluster"] = adata.obs["leiden"]

# Log cluster sizes
cluster_sizes = adata.obs['cluster'].value_counts().sort_index()
logger.info(f"Cluster sizes:\n{cluster_sizes.to_string()}")

# -----------------------------
# Plot UMAP with clusters and treatment metadata
# -----------------------------
logger.info("Generating UMAP plots colored by clusters and treatment")

# UMAP colored by cluster
sc.pl.umap(adata, color="cluster", save=f"_{FIG_NAME1}", show=False)

# UMAP colored by treatment/sample
if "treatment" in adata.obs.columns:
    FIG_NAME2 = "05_umap_treatment.png"
    sc.pl.umap(adata, color="treatment", save=f"_{FIG_NAME2}", show=False)
elif "sample" in adata.obs.columns:
    FIG_NAME2 = "05_umap_sample.png"
    sc.pl.umap(adata, color="sample", save=f"_{FIG_NAME2}", show=False)
else:
    logger.warning(
        "Metadata 'treatment' and 'sample' not found in adata.obs, skipping treatment/sample plot."
    )

# -----------------------------
# Marker Gene Identification
# -----------------------------
logger.info(f"Computing marker genes per cluster using {RANK_GENES_METHOD} method")
sc.tl.rank_genes_groups(adata, groupby="cluster", method=RANK_GENES_METHOD)

# Save marker gene results to CSV files
clusters = adata.obs["cluster"].cat.categories if hasattr(adata.obs["cluster"], "cat") else sorted(adata.obs["cluster"].unique())

summary_rows = []
for cluster in clusters:
    df = sc.get.rank_genes_groups_df(adata, group=cluster)
    df["cluster"] = cluster
    df = df[["cluster", "names", "logfoldchanges", "pvals_adj"]]
    df.columns = ["cluster", "gene", "logfoldchanges", "pvals_adj"]
    df.to_csv(MARKER_DIR / f"marker_genes_cluster_{cluster}.csv", index=False)
    logger.info(f"Saved marker genes for cluster {cluster} with {len(df)} genes")

    summary_rows.append(df.head(TOP_N_GENES))

# Save top marker genes summary
summary_df = pd.concat(summary_rows, ignore_index=True)
summary_df.to_csv(CSV_FILE, index=False)
logger.info(f"Saved top {TOP_N_GENES} marker genes summary for all clusters.")

# -----------------------------
# Generate HTML Report
# -----------------------------
figures = sorted(Path(sc.settings.figdir).glob("*05_*.png"))

with open(HTML_REPORT, "w") as f:
    f.write("<html><head><title>Step 5: Clustering & Marker Gene Analysis</title>\n")
    f.write("""
    <style>
        body { font-family: Arial, sans-serif; padding: 20px; }
        h1, h2 { color: #2c3e50; }
        img { border: 1px solid #ccc; padding: 10px; margin-bottom: 20px; width: 700px; }
        a { text-decoration: none; color: #2980b9; }
    </style>
    </head><body>
    """)
    f.write("<h1>Step 5 Report: Leiden Clustering and Marker Genes</h1>\n")

    f.write("<h2>UMAP Plots</h2>\n")
    for fig in figures:
        # Compute relative path from HTML_REPORT to the figure
        rel_path = os.path.relpath(fig, start=HTML_REPORT.parent)
        f.write(f'<img src="{rel_path}" alt="{fig.name}"><br>\n')

    f.write("<h2>Marker Genes per Cluster</h2>\n<ul>\n")
    for cluster in clusters:
        f.write(f'<li><a href="markers/marker_genes_cluster_{cluster}.csv">Marker genes for cluster {cluster}</a></li>\n')
    f.write("</ul>\n")

    f.write(f'<h2>Summary</h2>\n<p><a href="05_top_marker_genes_summary.csv">Top {TOP_N_GENES} marker genes summary (CSV)</a></p>\n')
    f.write("</body></html>")
logger.info(f"Generated HTML report at {HTML_REPORT}")

# -----------------------------
# Save the annotated data
# -----------------------------
adata.write(OUTPUT_FILE)
logger.info(f"Saved clustered AnnData object to {OUTPUT_FILE}")


print("✅ Step 5 complete: Clustering and marker gene analysis done and saved.")
