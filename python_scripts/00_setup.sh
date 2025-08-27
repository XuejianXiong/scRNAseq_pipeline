# -----------------------------
# Step 0: Setup paths
# -----------------------------

# setup required folders to store downloaded data, ouput files and figures

[ ! -d "data" ] && mkdir -p data

[ ! -d "data/cropseq" ] && mkdir -p data/cropseq

[ ! -d "data/cropseq/GSE149383" ] && mkdir -p data/cropseq/GSE149383
[ ! -d "data/cropseq/GSE149383/raw" ] && mkdir -p data/cropseq/GSE149383/raw
[ ! -d "data/cropseq/GSE149383/extracted" ] && mkdir -p data/cropseq/GSE149383/extracted

[ ! -d "data/cropseq/raw_fastq" ] && mkdir -p data/cropseq/raw_fastq

[ ! -d "data/genome" ] && mkdir -p data/genome

[ ! -d "data/retina" ] && mkdir -p data/retina
[ ! -d "data/retina/SRA559821" ] && mkdir -p data/retina/SRA559821
[ ! -d "data/retina/GSE137537" ] && mkdir -p data/retina/GSE137537

[ ! -d "results" ] && mkdir -p results
[ ! -d "results/cropseq" ] && mkdir -p results/cropseq
[ ! -d "results/retina" ] && mkdir -p results/retina

[ ! -d "results/retina/SRA559821" ] && mkdir -p results/retina/SRA559821
[ ! -d "results/retina/GSE137537" ] && mkdir -p results/retina/GSE137537

[ ! -d "figures" ] && mkdir -p figures
[ ! -d "figures/cropseq" ] && mkdir -p figures/cropseq
[ ! -d "figures/retina" ] && mkdir -p figures/retina

[ ! -d "figures/retina/SRA559821" ] && mkdir -p figures/retina/SRA559821
[ ! -d "figures/retina/GSE137537" ] && mkdir -p figures/retina/GSE137537
