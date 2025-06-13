# -----------------------------
# Step 2: Load and Merge Datasets with Metadata
# -----------------------------

library(Seurat)
library(Matrix)
library(dplyr)
library(SeuratObject)

# -----------------------------
# Parameters and Paths
# -----------------------------

data_dir <- "data/GSE149383"
output_rds <- "results/02_merged_data.rds"
output_meta <- "results/02_merged_metadata_R.csv"

if (!dir.exists("results")) dir.create("results")

sample_metadata <- list(
  "GSM3972651_PC9D0" = list(
    batch = "batch1",
    treatment = "DMSO",
    timepoint = "D0"
  ),
  "GSM3972652_PC9D3Erl" = list(
    batch = "batch1",
    treatment = "Erlotinib",
    timepoint = "D3"
  )
)

# -----------------------------
# Function to load a sample and add metadata
# -----------------------------

load_sample <- function(sample_id, meta) {
  cat(sprintf("📂 Loading: %s\n", sample_id))
  
  mat <- ReadMtx(
    mtx = file.path(data_dir, paste0(sample_id, "_matrix.mtx.gz")),
    features = file.path(data_dir, paste0(sample_id, "_features.tsv.gz")),
    cells = file.path(data_dir, paste0(sample_id, "_barcodes.tsv.gz"))
  )
  
  seurat_obj <- CreateSeuratObject(counts = mat, project = sample_id)
  
  # Add metadata columns
  for (key in names(meta)) {
    seurat_obj[[key]] <- meta[[key]]
  }
  seurat_obj[["sample"]] <- sample_id
  
  return(seurat_obj)
}

# -----------------------------
# Load all samples
# -----------------------------

seurat_list <- lapply(names(sample_metadata), function(id) {
  load_sample(id, sample_metadata[[id]])
})

# -----------------------------
# Merge Seurat objects
# -----------------------------

# Note: merge.data=TRUE merges metadata
seurat_merged <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), seurat_list)

cat(sprintf("✅ Merged data (multiple count layers): %d cells × %d genes\n", ncol(seurat_merged), nrow(seurat_merged)))

# -----------------------------
# Combine multiple "counts.*" layers into one "counts" layer
# -----------------------------

rna_assay <- seurat_merged[["RNA"]]
count_layers <- grep("^counts", SeuratObject::Layers(rna_assay), value = TRUE)

counts_list <- lapply(count_layers, function(layer) {
  SeuratObject::GetAssayData(rna_assay, layer = layer)
})

# Combine counts by columns (cells)
combined_counts <- do.call(cbind, counts_list)

# Create new RNA assay with combined counts
new_rna <- CreateAssayObject(counts = combined_counts)

# Replace old RNA assay with new combined assay
seurat_merged[["RNA"]] <- new_rna

# Confirm only one counts layer now
cat("Assay layers after combining counts:\n")
print(SeuratObject::Layers(seurat_merged[["RNA"]]))

cat(sprintf("✅ Combined counts into single 'counts' layer: %d cells × %d genes\n",
            ncol(seurat_merged), nrow(seurat_merged)))

# -----------------------------
# Save merged Seurat object and metadata
# -----------------------------

saveRDS(seurat_merged, file = output_rds)
write.csv(seurat_merged@meta.data, file = output_meta, row.names = TRUE)

cat("✅ Step 2 complete: Datasets loaded, merged, combined counts, and metadata saved.\n")
