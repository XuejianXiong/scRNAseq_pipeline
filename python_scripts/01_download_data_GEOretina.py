#!/usr/bin/env python3
"""
Prepare 10X-formatted input files for the retina dataset GSE137537.

This script:
1. Downloads the matrix, features, and metadata from GEO.
2. Formats features.tsv.gz (gene_id, gene_symbol, feature_type).
3. Extracts barcodes.tsv.gz from metadata.
4. Validates the matrix dimensions.
5. Renames files to standardized 10X prefixes.
"""

import gzip
import sys
import requests
import pandas as pd
from pathlib import Path
from scipy import io
from logzero import logger


# -----------------------------
# Configuration
# -----------------------------
DATA_NAME = "GSE137537"
FILES = {
    "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137537/suppl/GSE137537_counts.mtx.gz",
    "features": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137537/suppl/GSE137537_gene_names.txt.gz",
    "metadata": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137537/suppl/GSE137537_sample_annotations.tsv.gz",
}


# -----------------------------
# Helper functions
# -----------------------------
def download_file(url: str, dest: Path) -> None:
    """Download a file from URL to destination path."""
    logger.info(f"Downloading {url} → {dest}")
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        logger.error(f"Failed to download {url} (status {r.status_code})")
        sys.exit(1)
    with open(dest, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    logger.info(f"Downloaded {dest}")


def prepare_features(features_gz: Path, output_gz: Path) -> list[str]:
    """Create a 10X-style features.tsv.gz file."""
    logger.info("Preparing features.tsv.gz")
    with gzip.open(features_gz, "rt") as f:
        genes = [line.strip() for line in f]
    if not genes:
        logger.error("Gene list is empty.")
        sys.exit(1)

    rows = [f"{g}\t{g}\tGene Expression" for g in genes]
    with gzip.open(output_gz, "wt") as f:
        f.write("\n".join(rows))
    logger.info(f"Wrote {output_gz} ({len(genes)} genes)")
    return genes


def prepare_barcodes(meta_gz: Path, output_gz: Path) -> list[str]:
    """Create a barcodes.tsv.gz file using metadata."""
    logger.info("Preparing barcodes.tsv.gz")
    meta_df = pd.read_csv(meta_gz, sep="\t")
    if "Barcode" not in meta_df.columns:
        logger.error("Column 'Barcode' not found in metadata.")
        sys.exit(1)

    barcodes = meta_df["Barcode"].tolist()
    with gzip.open(output_gz, "wt") as f:
        f.write("\n".join(barcodes))
    logger.info(f"Wrote {output_gz} ({len(barcodes)} barcodes)")
    return barcodes


def verify_matrix(matrix_gz: Path, n_genes: int, n_barcodes: int) -> None:
    """Check that the matrix dimensions match features and barcodes."""
    logger.info("Verifying matrix.mtx.gz")
    try:
        with gzip.open(matrix_gz, "rb") as f:
            mat = io.mmread(f)
        if mat.shape[0] != n_genes or mat.shape[1] != n_barcodes:
            logger.error(
                f"Matrix dimensions mismatch: {mat.shape} vs "
                f"{n_genes} genes × {n_barcodes} barcodes"
            )
            sys.exit(1)
        logger.info(f"Matrix OK: {mat.shape[0]} genes × {mat.shape[1]} cells")
    except Exception as e:
        logger.error(f"Failed to read matrix: {e}")
        sys.exit(1)


# -----------------------------
# Main function
# -----------------------------
def main() -> None:
    outdir = Path("data") / "retina" / DATA_NAME
    outdir.mkdir(parents=True, exist_ok=True)

    matrix_gz = outdir / "GSE137537_counts.mtx.gz"
    features_gz = outdir / "GSE137537_gene_names.txt.gz"
    meta_gz = outdir / "GSE137537_sample_annotations.tsv.gz"

    # Step 1: Download GEO files
    if not matrix_gz.exists():
        download_file(FILES["matrix"], matrix_gz)
    if not features_gz.exists():
        download_file(FILES["features"], features_gz)
    if not meta_gz.exists():
        download_file(FILES["metadata"], meta_gz)

    # Step 2: Prepare features
    features_out = outdir / f"{DATA_NAME}_features.tsv.gz"
    genes = prepare_features(features_gz, features_out)

    # Step 3: Prepare barcodes
    barcodes_out = outdir / f"{DATA_NAME}_barcodes.tsv.gz"
    barcodes = prepare_barcodes(meta_gz, barcodes_out)

    # Step 4: Verify matrix
    verify_matrix(matrix_gz, len(genes), len(barcodes))

    # Step 5: Save matrix with standardized name
    matrix_out = outdir / f"{DATA_NAME}_matrix.mtx.gz"
    if not matrix_out.exists():
        matrix_out.write_bytes(matrix_gz.read_bytes())
    logger.info(f"Matrix saved to {matrix_out}")

    logger.info("✅ All 10X files prepared successfully.")


if __name__ == "__main__":
    main()
