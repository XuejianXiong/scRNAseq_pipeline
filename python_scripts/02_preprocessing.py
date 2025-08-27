# -----------------------------
# Step 2: Load datasets (single or multiple) based on SAMPLE_METADATA
# -----------------------------

import scanpy as sc
from logzero import logger
from pathlib import Path
import warnings

from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action='ignore', category=FutureWarning)

# -----------------------------
# User-Adjustable Parameters
# -----------------------------
#PROJ_NAME = "cropseq" 
PROJ_NAME = "retina"  

RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("02_log.txt", PROJ_NAME)
MERGED_DATA_FILE = RESULTS_DIR / "02_merged_data.h5ad"
METADATA_CSV = RESULTS_DIR / "02_merged_metadata.csv"

# Define metadata here; number of entries determines single vs multiple datasets
if PROJ_NAME == "cropseq":
    DATA_DIR = Path(f"data/{PROJ_NAME}/GSE149383")
    SAMPLE_METADATA = {
        "GSM3972651_PC9D0": {"batch": "batch1", "treatment": "DMSO", "timepoint": "D0"},
        "GSM3972652_PC9D3Erl": {"batch": "batch1", "treatment": "Erlotinib", "timepoint": "D3"}
    }
elif PROJ_NAME == "retina":
    DATA_DIR = Path(f"data/{PROJ_NAME}/SRA559821")
    SAMPLE_METADATA = {
        "SRA559821": {
            "tissue": "Retina", 
            "protocol": "C1_Fluidigm", 
            "species": "Homo_sapiens", 
            "instrument": "Illumina_HiSeq2500"
        }
    }

logger.info("Step 2 started: Load datasets based on SAMPLE_METADATA")


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
    logger.info(f"Loaded sample '{sample_id}' with {adata.n_obs} cells and {adata.n_vars} genes")

if len(SAMPLE_METADATA) > 1:
    # -----------------------------
    # Concatenate all samples into one AnnData object
    # -----------------------------
    adata_merged = adatas[0].concatenate(
        adatas[1:],
        batch_key="sample_id",
        batch_categories=list(SAMPLE_METADATA.keys())
    )
    logger.info(f"Merged {len(adatas)} samples: resulting data has {adata_merged.n_obs} cells × {adata_merged.n_vars} genes")

elif len(SAMPLE_METADATA) == 1:
    adata_merged = adatas[0]


# -----------------------------
# Save AnnData and metadata CSV
# -----------------------------
adatas[0].write(MERGED_DATA_FILE)
adatas[0].obs.to_csv(METADATA_CSV)
logger.info(f"Saved AnnData to '{MERGED_DATA_FILE}'")
logger.info(f"Saved metadata table to '{METADATA_CSV}'")
    

print("✅ Step 2 complete: Datasets loaded, merged, and metadata added.")