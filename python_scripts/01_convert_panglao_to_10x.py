#!/usr/bin/env python3
"""
Convert PanglaoDB tab-delimited matrix to 10X Genomics gzipped format.
Generates matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz.
"""

import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import os
from logzero import logger
import gzip

# === User parameters ===
sample = "SRA559821"
data_dir = f"data/retina/{sample}"  # Where to save 10X files
input_file = f"{data_dir}/{sample}.mat.gz"  # PanglaoDB matrix file

# === Step 0. Ensure output directory exists ===
os.makedirs(data_dir, exist_ok=True)

# === Step 1. Load the tab-delimited count matrix ===
logger.info(f"Loading PanglaoDB matrix: {input_file}")
df = pd.read_csv(input_file, sep="\t", compression="gzip", index_col=0)
logger.info(f"Matrix shape: {df.shape[0]} genes x {df.shape[1]} cells")

# === Step 2. Prepare 10X components ===
genes = df.index.astype(str)
barcodes = df.columns.astype(str)
count_matrix = sp.csc_matrix(df.values, dtype=int)

# === Step 3. Write features.tsv.gz ===
features_path = os.path.join(data_dir, f"{sample}_features.tsv.gz")
with gzip.open(features_path, "wt") as f:
    for gene in genes:
        f.write(f"{gene}\t{gene}\tGene Expression\n")
logger.info(f"Saved features.tsv.gz to {features_path}")

# === Step 4. Write barcodes.tsv.gz ===
barcodes_path = os.path.join(data_dir, f"{sample}_barcodes.tsv.gz")
with gzip.open(barcodes_path, "wt") as f:
    for bc in barcodes:
        f.write(f"{bc}\n")
logger.info(f"Saved barcodes.tsv.gz to {barcodes_path}")

# === Step 5. Write matrix.mtx.gz ===
matrix_path = os.path.join(data_dir, f"{sample}_matrix.mtx.gz")
with gzip.open(matrix_path, "wb") as f:
    mmwrite(f, count_matrix)
logger.info(f"Saved matrix.mtx.gz to {matrix_path}")

print("✅ Conversion complete. All 10X-style files are gzipped!")
