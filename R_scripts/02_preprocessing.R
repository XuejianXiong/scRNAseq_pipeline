# -----------------------------
# Step 2: Load and Merge Datasets with Metadata
# -----------------------------

# Clear environment, graphics, and console
rm(list = ls())            # Remove all variables from the workspace
graphics.off()             # Close all open graphics devices
cat("\014")                # Clear the console (works in RStudio)

library(Seurat)
library(Matrix)
library(dplyr)
library(SeuratObject)

source("R_scripts/pipeline_utils.R")

# -----------------------------
# Parameters and Paths
# -----------------------------

data_dir <- "data/GSE149383"
output_dir <- "results"
output_rds <- file.path(output_dir, "02_merged_data.rds")
output_meta <- file.path(output_dir, "02_merged_metadata_R.csv")
log_file <- file.path(output_dir, "02_log.txt")

if (!dir.exists(data_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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
# Functions
# -----------------------------
# Function to load a sample and add metadata
load_sample <- function(sample_id, meta) {
  cat(sprintf("📂 Loading: %s\n", sample_id))
  
  # Read the matrix, features, and barcodes files for this sample
  mat <- ReadMtx(
    mtx = file.path(data_dir, paste0(sample_id, "_matrix.mtx.gz")),     # Path to matrix file
    features = file.path(data_dir, paste0(sample_id, "_features.tsv.gz")), # Path to features file
    cells = file.path(data_dir, paste0(sample_id, "_barcodes.tsv.gz"))     # Path to barcodes file
  )
  
  # Create a Seurat object using the loaded matrix and assign the project name as sample_id
  seurat_obj <- CreateSeuratObject(counts = mat, project = sample_id)
  
  # Add metadata from the meta list as new columns in the Seurat object
  for (key in names(meta)) {
    seurat_obj[[key]] <- meta[[key]]
  }
  
  # Add an explicit 'sample' column to record the sample ID
  seurat_obj[["sample"]] <- sample_id
  
  return(seurat_obj)
}


log_msg("Step 2 started: Load and merge datasets with metadata")

# -----------------------------
# Load all samples into Seurat objects
# -----------------------------
# Iterate over each sample ID in sample_metadata and load its data + metadata
# Result: a list of Seurat objects, one for each sample
seurat_list <- lapply(names(sample_metadata), function(id) {
  load_sample(id, sample_metadata[[id]])
})

# -----------------------------
# Merge Seurat objects into a single Seurat object
# -----------------------------
# Use Reduce + merge to iteratively merge all Seurat objects in the list
# merge.data = TRUE ensures that metadata columns are combined properly
seurat_merged <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), seurat_list)
log_msg(sprintf("✅ Merged data (multiple count layers): %d cells × %d genes\n", ncol(seurat_merged), nrow(seurat_merged)))

# -----------------------------
# Combine multiple "counts.*" layers into one unified "counts" layer
# -----------------------------
# Extract the RNA assay
rna_assay <- seurat_merged[["RNA"]]

# Identify all count layers (layers named starting with "counts")
count_layers <- grep("^counts", SeuratObject::Layers(rna_assay), value = TRUE)

# Retrieve the matrix from each count layer
counts_list <- lapply(count_layers, function(layer) {
  SeuratObject::GetAssayData(rna_assay, layer = layer)
})

# Concatenate the matrices by columns (cells)
combined_counts <- do.call(cbind, counts_list)

# Create a new RNA assay with the combined counts matrix
new_rna <- CreateAssayObject(counts = combined_counts)

# Replace the old RNA assay with this new combined assay
seurat_merged[["RNA"]] <- new_rna

log_msg("Assay layers after combining counts:\n")
print(SeuratObject::Layers(seurat_merged[["RNA"]]))

log_msg(sprintf("✅ Combined counts into single 'counts' layer: %d cells × %d genes\n",
            ncol(seurat_merged), nrow(seurat_merged)))

# -----------------------------
# Save merged Seurat object and its metadata
# -----------------------------
saveRDS(seurat_merged, file = output_rds)
write.csv(seurat_merged@meta.data, file = output_meta, row.names = TRUE)


cat("✅ Step 2 complete: Datasets loaded, merged, combined counts, and metadata saved.\n")
