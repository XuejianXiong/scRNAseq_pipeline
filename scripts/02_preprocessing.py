import scanpy as sc
import pandas as pd
import os

# -----------------------------
# Step 2: Load and Merge Datasets with Metadata
# -----------------------------

# Define data directory and sample-specific metadata
data_dir = "data/GSE149383"

# Define sample metadata manually (can also be read from a CSV)
sample_metadata = {
    "GSM3972651_PC9D0": {
        "batch": "batch1",
        "treatment": "DMSO",
        "timepoint": "D0"
    },
    "GSM3972652_PC9D3Erl": {
        "batch": "batch1",
        "treatment": "Erlotinib",
        "timepoint": "D3"
    }
}

# Load each sample and assign metadata
adatas = []
for sample_id, meta in sample_metadata.items():
    adata = sc.read_10x_mtx(
        data_dir,
        var_names='gene_symbols',
        cache=True,
        prefix=sample_id + "_"
    )
    
    # Annotate metadata for each cell
    n_cells = adata.n_obs
    for key, value in meta.items():
        adata.obs[key] = value
    adata.obs["sample"] = sample_id

    adatas.append(adata)

# Concatenate all samples
adata_merged = adatas[0].concatenate(
    adatas[1:], 
    batch_key="sample_id", 
    batch_categories=list(sample_metadata.keys())
)

# Save merged data
os.makedirs("results", exist_ok=True)
adata_merged.write("results/02_merged_data.h5ad")

# Optional: Save metadata table
adata_merged.obs.to_csv("results/02_merged_metadata.csv")

print("✅ Step 2 complete: Datasets loaded, merged, and metadata added.")
