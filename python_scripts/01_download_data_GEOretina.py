#!/usr/bin/env python3

import gzip
import sys
import pandas as pd
import numpy as np
import requests
from pathlib import Path
from scipy import io
from logzero import logger

# -----------------------------
# User settings
# -----------------------------
DATA_NAME = "GSE137537"

OUTDIR = Path(f"data/retina/{DATA_NAME}")
OUTDIR.mkdir(parents=True, exist_ok=True)

# GEO file URLs
FILES = {
    "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137537/suppl/GSE137537_counts.mtx.gz",
    "features": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137537/suppl/GSE137537_gene_names.txt.gz",
    "metadata": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137537/suppl/GSE137537_sample_annotations.tsv.gz"
}

matrix_gz = OUTDIR / "GSE137537_counts.mtx.gz"
features_gz = OUTDIR / "GSE137537_gene_names.txt.gz"
meta_gz = OUTDIR / "GSE137537_sample_annotations.tsv.gz"

# -----------------------------
# Helper function: download a file
# -----------------------------
def download_file(url, dest):
    logger.info(f"Downloading {url} -> {dest}")
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        logger.error(f"Failed to download {url}")
        sys.exit(1)
    with open(dest, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    logger.info(f"Downloaded {dest}")

# -----------------------------
# Step 1: Download GEO files
# -----------------------------

if not matrix_gz.exists():
    download_file(FILES["matrix"], matrix_gz)
if not features_gz.exists():
    download_file(FILES["features"], features_gz)
if not meta_gz.exists():
    download_file(FILES["metadata"], meta_gz)

# -----------------------------
# Step 2: Prepare features.tsv.gz
# -----------------------------
logger.info("Preparing features.tsv.gz")
with gzip.open(features_gz, "rt") as f:
    genes = [line.strip() for line in f]
if not genes:
    logger.error("Gene list is empty.")
    sys.exit(1)

# 10X expects 3 columns: gene_id, gene_symbol, feature_type
feature_rows = [f"{g}\t{g}\tGene Expression" for g in genes]
features_out = OUTDIR / f"{DATA_NAME}_features.tsv.gz"
with gzip.open(features_out, "wt") as f:
    f.write("\n".join(feature_rows))
logger.info(f"Wrote {features_out} ({len(genes)} genes)")

# -----------------------------
# Step 3: Prepare barcodes.tsv.gz
# -----------------------------
logger.info("Preparing barcodes.tsv.gz")
meta_df = pd.read_csv(meta_gz, sep="\t")
if "Barcode" not in meta_df.columns:
    logger.error("Column 'Barcode' not found in metadata.")
    sys.exit(1)

barcodes = meta_df["Barcode"].tolist()
barcodes_out = OUTDIR / f"{DATA_NAME}_barcodes.tsv.gz"
with gzip.open(barcodes_out, "wt") as f:
    f.write("\n".join(barcodes))
logger.info(f"Wrote {barcodes_out} ({len(barcodes)} barcodes)")

# -----------------------------
# Step 4: Verify matrix.mtx.gz
# -----------------------------
logger.info("Verifying matrix.mtx.gz")
try:
    with gzip.open(matrix_gz, "rb") as f:
        mat = io.mmread(f)
    if mat.shape[1] != len(barcodes) or mat.shape[0] != len(genes):
        logger.error("Matrix dimensions do not match gene/barcode lists.")
        sys.exit(1)
    logger.info(f"Matrix OK: {mat.shape[0]} genes × {mat.shape[1]} cells")
except Exception as e:
    logger.error(f"Failed to read matrix: {e}")
    sys.exit(1)

# -----------------------------
# Step 5: Save matrix with proper prefix (already gzipped)
# -----------------------------
matrix_out = OUTDIR / f"{DATA_NAME}_matrix.mtx.gz"
if not matrix_out.exists():
    # just rename/copy the downloaded matrix
    matrix_out.write_bytes(matrix_gz.read_bytes())
logger.info(f"Matrix saved to {matrix_out}")

logger.info("✅ All 10X files prepared successfully.")
