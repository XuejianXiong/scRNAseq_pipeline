# -----------------------------
# Step 2: Load and Merge Datasets with Metadata
# -----------------------------

import scanpy as sc
import pandas as pd
import logging
from pathlib import Path
import os
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# User-Adjustable Parameters
# -----------------------------

DATA_DIR = Path("data/GSE149383")
RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("02_log.txt")
MERGED_DATA_FILE = RESULTS_DIR / "02_merged_data.h5ad"
METADATA_CSV = RESULTS_DIR / "02_merged_metadata.csv"

# Sample metadata dictionary (manually defined)
SAMPLE_METADATA = {
    "GSM3972651_PC9D0": {"batch": "batch1", "treatment": "DMSO", "timepoint": "D0"},
    "GSM3972652_PC9D3Erl": {"batch": "batch1", "treatment": "Erlotinib", "timepoint": "D3"}
}
# -----------------------------


logging.info("Step 2 started: Load and merge datasets with metadata")

# -----------------------------
# Load each sample, annotate metadata, and collect AnnData objects
# -----------------------------
adatas = []
for sample_id, meta in SAMPLE_METADATA.items():
    adata = sc.read_10x_mtx(
        DATA_DIR,
        var_names='gene_symbols',
        cache=True,
        prefix=f"{sample_id}_"
    )
    for key, value in meta.items():
        adata.obs[key] = value
    adata.obs["sample"] = sample_id

    adatas.append(adata)
    logging.info(f"Loaded sample '{sample_id}' with {adata.n_obs} cells and {adata.n_vars} genes")

# -----------------------------
# Concatenate all samples into one AnnData object
# -----------------------------
adata_merged = adatas[0].concatenate(
    adatas[1:],
    batch_key="sample_id",
    batch_categories=list(SAMPLE_METADATA.keys())
)
logging.info(f"Merged {len(adatas)} samples: resulting data has {adata_merged.n_obs} cells × {adata_merged.n_vars} genes")

# -----------------------------
# Save merged AnnData and metadata CSV
# -----------------------------
adata_merged.write(MERGED_DATA_FILE)
adata_merged.obs.to_csv(METADATA_CSV)
logging.info(f"Saved merged AnnData to '{MERGED_DATA_FILE}'")
logging.info(f"Saved merged metadata table to '{METADATA_CSV}'")


print("✅ Step 2 complete: Datasets loaded, merged, and metadata added.")
