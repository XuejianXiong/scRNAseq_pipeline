# -----------------------------
# Step 0: Setup paths
# -----------------------------

# setup required folders to store downloaded data, ouput files and figures

[ ! -d "data" ] && mkdir -p data

[ ! -d "data/GSE149383" ] && mkdir -p data/GSE149383
[ ! -d "data/GSE149383/raw" ] && mkdir -p data/GSE149383/raw
[ ! -d "data/GSE149383/extracted" ] && mkdir -p data/GSE149383/extracted

[ ! -d "data/genome" ] && mkdir -p data/genome
[ ! -d "data/raw_fastq" ] && mkdir -p data/raw_fastq


[ ! -d "results" ] && mkdir -p results

[ ! -d "figures" ] && mkdir -p figures