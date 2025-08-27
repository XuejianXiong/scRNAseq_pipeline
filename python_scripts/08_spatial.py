import scanpy as sc
import squidpy as sq
from pathlib import Path
import os
import logging


# ----------------------------
# Config
# ----------------------------
# Input: processed retina scRNA-seq AnnData
scRNA_adata_path = "results/retina/adata_retina_scRNA.h5ad"

# Output directory
outdir = Path("results/retina/spatial")
outdir.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Step 1: Load spatial transcriptomics dataset
# Using Squidpy demo Visium dataset (replace with 10x human retina if desired)
adata_sp = sq.datasets.visium_hne_adata()  # small demo dataset

# If you downloaded real 10x retina Visium dataset:
# adata_sp = sc.read_visium(path="data/retina/spatial/")

# ----------------------------
# Step 2: Preprocessing
# ----------------------------
sc.pp.calculate_qc_metrics(adata_sp, inplace=True)
sc.pp.filter_genes(adata_sp, min_counts=3)
sc.pp.normalize_total(adata_sp, target_sum=1e4)
sc.pp.log1p(adata_sp)

# Compute neighbors & UMAP
sc.pp.highly_variable_genes(adata_sp, flavor="seurat_v3", n_top_genes=2000)
sc.pp.pca(adata_sp)
sc.pp.neighbors(adata_sp)
sc.tl.umap(adata_sp)

# ----------------------------
# Step 3: Load scRNA clusters
# ----------------------------
scRNA_adata = sc.read_h5ad(scRNA_adata_path)
# Assume cluster labels are stored in scRNA_adata.obs['leiden']

# ----------------------------
# Step 4: Project cluster signatures
# ----------------------------
# Compute top marker genes per cluster
marker_dict = {}
for cluster in scRNA_adata.obs['leiden'].cat.categories:
    df = sc.tl.rank_genes_groups(scRNA_adata, groupby='leiden', method='t-test', n_genes=50, pts=True, copy=True)
    # Use top 30 genes per cluster
    top_genes = df['names'][cluster][:30]
    marker_dict[cluster] = list(top_genes)

# Score gene programs on spatial data
for cluster, genes in marker_dict.items():
    # Keep only genes present in spatial dataset
    genes_filtered = [g for g in genes if g in adata_sp.var_names]
    if not genes_filtered:
        print(f"Warning: No overlapping genes for cluster {cluster}")
        continue
    sc.tl.score_genes(adata_sp, gene_list=genes_filtered, score_name=f"cluster_{cluster}_score")
    # Plot spatial distribution
    sc.pl.spatial(
        adata_sp,
        color=f"cluster_{cluster}_score",
        spot_size=1.5,
        cmap='viridis',
        title=f"Cluster {cluster} program",
        save=f"_cluster_{cluster}_score.png"
    )

# ----------------------------
# Step 5: Save processed spatial object
# ----------------------------
adata_sp.write_h5ad(outdir / "adata_retina_spatial_scored.h5ad")
print("Spatial transcriptomics mapping completed. Figures saved to", outdir)
