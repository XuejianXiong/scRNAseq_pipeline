import scanpy as sc
import pandas as pd
import logging
from pathlib import Path
import os

# -----------------------------
# Step 5: Clustering and Marker Gene Analysis
# -----------------------------

# -----------------------------
# User-adjustable parameters
# -----------------------------
LEIDEN_RESOLUTION = 0.5        # Adjust clustering granularity (higher = more clusters)
N_NEIGHBORS = 15               # Number of neighbors to use for graph construction (adjust based on data)
N_PCS = 40                     # Number of PCs to use
RANK_GENES_METHOD = "wilcoxon"  # Marker gene ranking method
TOP_N_GENES = 10               # Number of top marker genes per cluster to save in summary

# -----------------------------
# Setup
# -----------------------------
sc.settings.verbosity = 2
sc.settings.autoshow = False
sc.settings.figdir = "figures"
Path(sc.settings.figdir).mkdir(parents=True, exist_ok=True)

os.makedirs("results", exist_ok=True)
log_file = "results/05_log.txt"
logging.basicConfig(filename=log_file, level=logging.INFO, format="%(asctime)s %(message)s")
logging.info("Step 5 started: Clustering and Marker Gene Analysis")

logging.info(f"Parameters - LEIDEN_RESOLUTION: {LEIDEN_RESOLUTION}, N_NEIGHBORS: {N_NEIGHBORS}, "
             f"N_PCS: {N_PCS}, RANK_GENES_METHOD: {RANK_GENES_METHOD}, TOP_N_GENES: {TOP_N_GENES}")

# Load normalized and reduced data from Step 4
adata = sc.read("results/04_norm_dr_data.h5ad")
#print(adata.raw)  # Should not be None
#print(adata.raw.X[:5, :5])  # Check some values (should look log-transformed)
#print(adata.X[:5, :5])      # Check current data matrix values for comparison
logging.info(f"Loaded normalized data: {adata.n_obs} cells, {adata.n_vars} genes")

# -----------------------------
# Recompute neighbors if needed (adjust parameters)
# -----------------------------
logging.info(f"Computing neighbors with n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}")
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)

# -----------------------------
# Clustering with Leiden algorithm
# -----------------------------
logging.info(f"Running Leiden clustering with resolution={LEIDEN_RESOLUTION}")
sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION)

# Save clustering labels in obs
adata.obs["cluster"] = adata.obs["leiden"]

# Log cluster sizes
cluster_sizes = adata.obs['cluster'].value_counts().sort_index()
logging.info(f"Cluster sizes:\n{cluster_sizes.to_string()}")

# -----------------------------
# Plot UMAP with clusters and treatment metadata
# -----------------------------
logging.info("Generating UMAP plots colored by clusters and treatment")

# UMAP colored by cluster
sc.pl.umap(adata, color="cluster", save="_05_umap_clusters.png", show=False)

# UMAP colored by treatment (if present)
if "treatment" in adata.obs.columns:
    sc.pl.umap(adata, color="treatment", save="_05_umap_treatment.png", show=False)
else:
    logging.warning("Metadata 'treatment' not found in adata.obs, skipping treatment plot.")

# -----------------------------
# Marker Gene Identification
# -----------------------------
logging.info(f"Computing marker genes per cluster using {RANK_GENES_METHOD} method")
sc.tl.rank_genes_groups(adata, groupby="cluster", method=RANK_GENES_METHOD)

# Save marker gene results to CSV files
marker_dir = Path("results/markers")
marker_dir.mkdir(exist_ok=True)

clusters = adata.obs["cluster"].cat.categories if hasattr(adata.obs["cluster"], "cat") else sorted(adata.obs["cluster"].unique())
summary_rows = []

for cluster in clusters:
    df = sc.get.rank_genes_groups_df(adata, group=cluster)
    df["cluster"] = cluster
    df = df[["cluster", "names", "logfoldchanges", "pvals_adj"]]
    df.columns = ["cluster", "gene", "logfoldchanges", "pvals_adj"]
    df.to_csv(marker_dir / f"marker_genes_cluster_{cluster}.csv", index=False)
    logging.info(f"Saved marker genes for cluster {cluster} with {len(df)} genes")

    summary_rows.append(df.head(TOP_N_GENES))

# Save top marker genes summary
summary_df = pd.concat(summary_rows, ignore_index=True)
summary_df.to_csv("results/05_top_marker_genes_summary.csv", index=False)
logging.info(f"Saved top {TOP_N_GENES} marker genes summary for all clusters.")

# -----------------------------
# Generate HTML Report
# -----------------------------
report_path = Path("results/05_report.html")
figures = sorted(Path(sc.settings.figdir).glob("*05_*.png"))

with open(report_path, "w") as f:
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
        f.write(f'<img src="../{fig}" alt="{fig.name}"><br>\n')

    f.write("<h2>Marker Genes per Cluster</h2>\n<ul>\n")
    for cluster in clusters:
        f.write(f'<li><a href="../results/markers/marker_genes_cluster_{cluster}.csv">Marker genes for cluster {cluster}</a></li>\n')
    f.write("</ul>\n")

    f.write(f'<h2>Summary</h2>\n<p><a href="../results/05_top_marker_genes_summary.csv">Top {TOP_N_GENES} marker genes summary (CSV)</a></p>\n')
    f.write("</body></html>")

logging.info(f"Generated HTML report at {report_path}")

# -----------------------------
# Save the annotated data
# -----------------------------
adata.write("results/05_clustered_data.h5ad")
logging.info("Saved clustered AnnData object to results/05_clustered_data.h5ad")

print("✅ Step 5 complete: Clustering and marker gene analysis done and saved.")
