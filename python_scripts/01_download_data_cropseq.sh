#!/bin/bash
set -euo pipefail

# -----------------------------
# Directories
# -----------------------------
CROPSEQ_RAW_DIR="data/cropseq/GSE149383/raw"
CROPSEQ_EXTRACT_DIR="data/cropseq/GSE149383/extracted"
FASTQ_DIR="data/cropseq/raw_fastq"
GENOME_DIR="data/genome"

mkdir -p "$CROPSEQ_RAW_DIR" "$CROPSEQ_EXTRACT_DIR" "$FASTQ_DIR" "$GENOME_DIR"

# -----------------------------
# Download and extract GEO data
# -----------------------------
cd "$CROPSEQ_RAW_DIR"

GSE_TAR="GSE149383_RAW.tar"

if [ ! -f "$GSE_TAR" ]; then
    echo "Downloading GSE149383_RAW.tar..."
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149383/suppl/$GSE_TAR
else
    echo "$GSE_TAR already exists, skipping download."
fi

# Extract main tar archive if extraction directory is empty
if [ -z "$(ls -A ../extracted 2>/dev/null)" ]; then
    echo "Extracting $GSE_TAR..."
    mkdir -p ../extracted
    tar -xvf "$GSE_TAR" -C ../extracted
else
    echo "GSE149383 already extracted, skipping extraction."
fi

# Process inner tar.gz files
cd ../extracted
for archive in GSM397265*_filtered_feature_bc_matrices.tar.gz; do
    prefix="${archive%%_filtered_feature_bc_matrices.tar.gz}"
    TEMP_DIR="temp_extract"

    if [ ! -f "../${prefix}_barcodes.tsv" ]; then
        echo "Processing $archive..."
        mkdir -p "$TEMP_DIR"
        tar -xzf "$archive" -C "$TEMP_DIR"
        find "$TEMP_DIR" -type f | while read filepath; do
            filename=$(basename "$filepath")
            cp "$filepath" "../${prefix}_${filename}"
        done
        rm -rf "$TEMP_DIR"
    else
        echo "$prefix files already exist, skipping."
    fi
done

cd ../../../../  # back to project root

# -----------------------------
# Download raw FASTQ files
# -----------------------------
cd "$FASTQ_DIR"
FASTQ_TAR="pbmc_1k_v3_fastqs.tar"

if [ ! -f "$FASTQ_TAR" ]; then
    echo "Downloading PBMC 1k v3 FASTQ files..."
    wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/$FASTQ_TAR
else
    echo "$FASTQ_TAR already exists, skipping download."
fi

# Extract if FASTQ files not present
if [ ! -d "pbmc_1k_v3" ]; then
    echo "Extracting FASTQ files..."
    tar -xvf "$FASTQ_TAR"
else
    echo "FASTQ files already extracted, skipping."
fi

cd ../../../  # back to project root

# -----------------------------
# Download hg38 reference
# -----------------------------
cd "$GENOME_DIR"

FASTA="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="Homo_sapiens.GRCh38.109.gtf"

if [ ! -f "$FASTA" ]; then
    echo "Downloading hg38 FASTA..."
    wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/$FASTA.gz
    gunzip $FASTA.gz
else
    echo "$FASTA already exists, skipping."
fi

if [ ! -f "$GTF" ]; then
    echo "Downloading hg38 GTF..."
    wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/$GTF.gz
    gunzip $GTF.gz
else
    echo "$GTF already exists, skipping."
fi

cd ../../  # back to project root

echo "All downloads and extractions complete."
