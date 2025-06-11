cd data/GSE149383/raw

# Download raw tar archive (alternatively use wget or curl)
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149383/suppl/GSE149383_RAW.tar

#tar -xvf GSE149383_RAW.tar -C ../extracted

cd ../extracted

# List of your tar.gz files
for archive in GSM397265*_filtered_feature_bc_matrices.tar.gz; do
    # Extract prefix from filename (remove .tar.gz)
    prefix="${archive%%_filtered_feature_bc_matrices.tar.gz}"

    # Create a temp directory
    mkdir -p temp_extract
    tar -xzf "$archive" -C temp_extract

    # Find all files (not directories) in extracted folder
    find temp_extract -type f | while read filepath; do
        filename=$(basename "$filepath")
        cp "$filepath" "../${prefix}_${filename}"
    done

    # Cleanup
    rm -rf temp_extract
done
