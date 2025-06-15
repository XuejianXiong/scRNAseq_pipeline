# -----------------------------
# Step 4: Normalization and Dimensionality Reduction
# -----------------------------

# Clear environment, graphics, and console
rm(list = ls())            # Remove all variables from the workspace
graphics.off()             # Close all open graphics devices
cat("\014")                # Clear the console (works in RStudio)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# -----------------------------
# Parameters and Paths
# -----------------------------
output_dir <- "results"
plot_dir <- "figures"
input_file <- file.path(output_dir, "03_filtered_data.rds")
output_file <- file.path(output_dir, "04_norm_dr_data.rds")
hvg_file <- file.path(plot_dir, "04_hvg_plot.png")
pca_file <- file.path(plot_dir, "04_pca_treatment.png")
umap_file <- file.path(plot_dir, "04_umap_treatment.png")
report_file <- file.path(output_dir, "04_report.html")
log_file <- file.path(output_dir, "04_log.txt")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

norm_params <- list(
  n_top_hvg = 2000,
  scale_max = 10,
  n_pcs = 40,
  n_neighbors = 15
)


# -----------------------------
# Functions
# -----------------------------
# Function to set the log message
log_msg <- function(msg) {
  
  # Generate current timestamp in YYYY-MM-DD HH:MM:SS format
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Construct log line with timestamp prefix
  line <- sprintf("[%s] %s", timestamp, msg)
  
  # Print log message to console
  cat(line, "\n")
  
  # Append log message to the specified log file
  # Assumes 'log_file' is a global variable set elsewhere in the script
  write(line, file = log_file, append = TRUE)
}


log_msg("Step 4 started: Normalization and Dimensionality Reduction")

# -----------------------------
# Load Data
# -----------------------------
seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded filtered data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# Normalization and HVG Selection
# -----------------------------
# Normalize gene expression using LogNormalize method: 
# each cell's counts are scaled to a total of 10,000, then log-transformed
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
log_msg("Normalized total counts per cell to 10,000")

# Identify highly variable genes (HVGs) 
# using variance-stabilizing transformation (VST) method
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = norm_params$n_top_hvg)
log_msg(sprintf("Selected top %d highly variable genes", norm_params$n_top_hvg))

# -----------------------------
# HVG Plot
# -----------------------------
# Generate scatter plot showing variance of genes, highlight top HVGs
hvg_plot <- VariableFeaturePlot(seurat_obj)

# Label the top 10 variable genes on the plot for reference
hvg_plot <- LabelPoints(
  plot = hvg_plot,
  points = head(VariableFeatures(seurat_obj), 10),
  repel = TRUE
)

# Save HVG plot as an image file
ggsave(hvg_file, hvg_plot, width = 7, height = 5)
log_msg("Saved HVG plot")

# -----------------------------
# Scaling, PCA, and UMAP
# -----------------------------
# Scale data (center and scale HVG features); 
# clip extreme values for stability if specified
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), clip.max = norm_params$scale_max)
log_msg(sprintf("Scaled data with max clip value of %d", norm_params$scale_max))

# Perform PCA on scaled HVG data to reduce dimensionality
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = norm_params$n_pcs)
log_msg(sprintf("Performed PCA using %d PCs", norm_params$n_pcs))

# -----------------------------
# PCA Plot
# -----------------------------
# Plot PCA results, color points by treatment group
p1 <- DimPlot(seurat_obj, reduction = "pca", group.by = "treatment")

# Save PCA plot to file
ggsave(pca_file, p1, width = 7, height = 5)
log_msg("Saved PCA plot colored by treatment")

# -----------------------------
# Neighbors and UMAP
# -----------------------------
# Compute k-nearest neighbors in PCA space for clustering/UMAP
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:norm_params$n_pcs, k.param = norm_params$n_neighbors)

# Run UMAP for non-linear dimensionality reduction
seurat_obj <- RunUMAP(seurat_obj, dims = 1:norm_params$n_pcs)
log_msg(sprintf("Computed neighbors (k=%d) and UMAP", norm_params$n_neighbors))

# -----------------------------
# UMAP Plot
# -----------------------------
# Plot UMAP embedding, color by treatment group
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "treatment")

# Save UMAP plot to file
ggsave(umap_file, p2, width = 7, height = 5)
log_msg("Saved UMAP plot colored by treatment")

# -----------------------------
# Save Processed Data
# -----------------------------
# Save processed Seurat object (normalized, HVGs, PCA, UMAP) for downstream steps
saveRDS(seurat_obj, output_file)
log_msg(sprintf("Saved normalized and reduced data to %s", output_file))


# -----------------------------
# HTML Report
# -----------------------------
# Initialize HTML content as a character vector (each element = one line)
html <- c(
  # HTML header with document title
  "<html><head><title>Normalization & Dimensionality Reduction Report</title></head><body>",
  
  # Main heading for the report
  "<h1>Step 4 Report: Normalization, HVG, PCA, UMAP</h1>",
  
  # Insert timestamp for when the report was generated
  sprintf("<p><b>Date:</b> %s</p>", Sys.time()),
  
  # Document the key parameters used in this analysis step for transparency and reproducibility
  sprintf("<p><b>Parameters:</b> n_top_hvg = %d, scale_max = %d, n_pcs = %d, n_neighbors = %d</p>",
          norm_params$n_top_hvg, norm_params$scale_max, norm_params$n_pcs, norm_params$n_neighbors),
  
  # Section for Highly Variable Genes with linked HVG plot
  sprintf('<h2>Highly Variable Genes</h2><img src="../%s/04_hvg_plot.png" width="600"/><br>', plot_dir),
  
  # Section for PCA results with linked PCA plot (colored by treatment)
  sprintf('<h2>PCA</h2><img src="../%s/04_pca_treatment.png" width="600"/><br>', plot_dir),
  
  # Section for UMAP embedding with linked UMAP plot (colored by treatment)
  sprintf('<h2>UMAP</h2><img src="../%s/04_umap_treatment.png" width="600"/><br>', plot_dir),
  
  # Close HTML body and document
  "</body></html>"
)

# Write the assembled HTML lines to file
writeLines(html, report_file)
log_msg(sprintf("Generated HTML report at: %s", report_file))


cat("✅ Step 4 complete: Normalization, HVG selection, PCA & UMAP done and saved.\n")
