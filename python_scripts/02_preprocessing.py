#!/usr/bin/env python3
"""
Step 2: Load single-cell datasets (single or multiple samples)
and merge into a single AnnData object with annotated metadata.
"""

import json
import warnings
from pathlib import Path
import argparse

import scanpy as sc
from logzero import logger
from pipeline_utils import setup_dirs_logs

warnings.simplefilter(action="ignore", category=FutureWarning)


# -----------------------------
# Helper Functions
# -----------------------------
def load_metadata_config(metadata_file: Path, proj_name: str):
    """
    Load dataset directory and sample metadata from a JSON config file.
    """
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata config file not found: {metadata_file}")

    with open(metadata_file, "r") as f:
        config = json.load(f)

    if proj_name not in config:
        raise ValueError(f"Project '{proj_name}' not found in {metadata_file}")

    data_dir = Path(config[proj_name]["data_dir"])
    sample_metadata = config[proj_name]["samples"]

    if not data_dir.exists():
        logger.warning(f"Data directory does not exist yet: {data_dir}")

    return data_dir, sample_metadata


def load_and_annotate_samples(data_dir: Path, sample_metadata: dict):
    """
    Load 10X datasets for each sample, annotate metadata,
    and return a list of AnnData objects.
    """
    adatas = []
    for sample_id, meta in sample_metadata.items():
        adata = sc.read_10x_mtx(
            data_dir,
            var_names="gene_symbols",
            cache=True,
            prefix=f"{sample_id}_",
        )
        for key, value in meta.items():
            adata.obs[key] = value
        adata.obs["sample"] = sample_id

        logger.info(
            f"Loaded sample '{sample_id}' with {adata.n_obs} cells × {adata.n_vars} genes"
        )
        adatas.append(adata)
    return adatas


def merge_adatas(adatas: list, sample_metadata: dict):
    """
    Merge a list of AnnData objects into a single AnnData if multiple samples exist.
    """
    if len(adatas) > 1:
        adata_merged = adatas[0].concatenate(
            adatas[1:],
            batch_key="sample_id",
            batch_categories=list(sample_metadata.keys()),
        )
        logger.info(
            f"Merged {len(adatas)} samples: resulting data has "
            f"{adata_merged.n_obs} cells × {adata_merged.n_vars} genes"
        )
    elif len(adatas) == 1:
        adata_merged = adatas[0]
    else:
        raise ValueError("No samples were loaded.")
    return adata_merged


def save_outputs(adata, results_dir: Path):
    """
    Save merged AnnData object and its metadata table.
    """
    merged_file = results_dir / "02_merged_data.h5ad"
    metadata_csv = results_dir / "02_merged_metadata.csv"

    adata.write(merged_file)
    adata.obs.to_csv(metadata_csv)

    logger.info(f"Saved AnnData to '{merged_file}'")
    logger.info(f"Saved metadata table to '{metadata_csv}'")


# -----------------------------
# Main Pipeline
# -----------------------------
def main():

    # Parameters:

    # PROJ_NAME = "cropseq" or "retina_SRA559821" or "retina_GSE137537"
    parser = argparse.ArgumentParser(
        description="Step 2: Load and merge scRNA-seq datasets"
    )
    parser.add_argument("project_name", help="Project name as defined in metadata.json")
    args = parser.parse_args()

    proj_name = args.project_name

    METADATA_FILE = Path("data/metadata_config.json")

    # Set up directories and logging
    RESULTS_DIR, FIGURE_DIR, LOG_FILE = setup_dirs_logs("02_log.txt", proj_name)
    logger.info("Step 2 started: Load datasets based on SAMPLE_METADATA")

    # Load sample metadata
    data_dir, sample_metadata = load_metadata_config(METADATA_FILE, proj_name)

    # Load and annotate samples
    adatas = load_and_annotate_samples(data_dir, sample_metadata)

    # Merge datasets if needed
    adata_merged = merge_adatas(adatas, sample_metadata)

    # Save outputs
    save_outputs(adata_merged, RESULTS_DIR)

    logger.info("✅ Step 2 complete: Datasets loaded, merged, and metadata added.")


if __name__ == "__main__":
    main()
