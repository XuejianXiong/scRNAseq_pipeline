
: << 'SKIP_BLOCK'
SKIP_BLOCK
# -----------------------------
# Download data from GEO
# -----------------------------
# Move to the directory for storing raw GSE149383 files
cd data/GSE149383/raw

# Download the raw tar archive containing supplementary data files
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149383/suppl/GSE149383_RAW.tar

# Extract the main tar archive into the ../extracted directory
tar -xvf GSE149383_RAW.tar -C ../extracted

# Move into the extracted directory to process inner tar.gz files
cd ../extracted

# Loop over each filtered_feature_bc_matrices tar.gz archive
for archive in GSM397265*_filtered_feature_bc_matrices.tar.gz; do
    # Extract prefix for naming outputs (remove _filtered_feature_bc_matrices.tar.gz)
    prefix="${archive%%_filtered_feature_bc_matrices.tar.gz}"

    # Create a temporary directory for extraction
    mkdir -p temp_extract
    
    # Extract the contents of the current tar.gz file into the temp directory
    tar -xzf "$archive" -C temp_extract

    # Copy all extracted files (not directories) to the parent directory
    # Prefix the filenames with the sample identifier for clarity
    find temp_extract -type f | while read filepath; do
        filename=$(basename "$filepath")
        cp "$filepath" "../${prefix}_${filename}"
    done

    # Clean up the temporary directory after processing each archive
    rm -rf temp_extract
done

# Return to the root directory of your project
cd ../../../


# -----------------------------
# Download raw FASTQ files
# -----------------------------
# Move to the directory where raw FASTQ files will be stored
cd data/raw_fastq

# Downloading 10X PBMC 1k v3 dataset
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar

# Extracting FASTQ files
tar -xvf pbmc_1k_v3_fastqs.tar

# Return to the project root directory
cd ../../


# -----------------------------
# Download hg38 reference for STARsolo
# -----------------------------
cd data/genome

# ownloading hg38 reference fasta
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Downloading hg38 GTF annotation
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Unzipping reference files
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz

# Return to the project root directory
cd ../../


echo "Download complete."
