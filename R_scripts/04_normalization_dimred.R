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
# Logging Utility
# -----------------------------

log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n")
  write(line, file = log_file, append = TRUE)
}

# -----------------------------
# Load Data
# -----------------------------

seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded filtered data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# Normalization and HVG Selection
# -----------------------------

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
log_msg("Normalized total counts per cell to 10,000")

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = norm_params$n_top_hvg)
log_msg(sprintf("Selected top %d highly variable genes", norm_params$n_top_hvg))

# Save HVG plot
hvg_plot <- VariableFeaturePlot(seurat_obj)
hvg_plot <- LabelPoints(
  plot = hvg_plot,
  points = head(VariableFeatures(seurat_obj), 10),
  repel = TRUE
)
ggsave(hvg_file, hvg_plot, width = 7, height = 5)
log_msg("Saved HVG plot")

# -----------------------------
# Scaling, PCA, and UMAP
# -----------------------------

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), clip.max = norm_params$scale_max)
log_msg(sprintf("Scaled data with max clip value of %d", norm_params$scale_max))

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = norm_params$n_pcs)
log_msg(sprintf("Performed PCA using %d PCs", norm_params$n_pcs))

# PCA plot
p1 <- DimPlot(seurat_obj, reduction = "pca", group.by = "treatment")
ggsave(pca_file, p1, width = 7, height = 5)
log_msg("Saved PCA plot colored by treatment")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:norm_params$n_pcs, k.param = norm_params$n_neighbors)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:norm_params$n_pcs)
log_msg(sprintf("Computed neighbors (k=%d) and UMAP", norm_params$n_neighbors))

# UMAP plot
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "treatment")
ggsave(umap_file, p2, width = 7, height = 5)
log_msg("Saved UMAP plot colored by treatment")

# -----------------------------
# Save Processed Data
# -----------------------------

saveRDS(seurat_obj, output_file)
log_msg(sprintf("Saved normalized and reduced data to %s", output_file))

# -----------------------------
# Simple HTML Report
# -----------------------------

html <- c(
  "<html><head><title>Normalization & Dimensionality Reduction Report</title></head><body>",
  "<h1>Step 4 Report: Normalization, HVG, PCA, UMAP</h1>",
  sprintf("<p><b>Date:</b> %s</p>", Sys.time()),
  sprintf("<p><b>Parameters:</b> n_top_hvg = %d, scale_max = %d, n_pcs = %d, n_neighbors = %d</p>",
          norm_params$n_top_hvg, norm_params$scale_max, norm_params$n_pcs, norm_params$n_neighbors),
  sprintf('<h2>Highly Variable Genes</h2><img src="../%s/04_hvg_plot.png" width="600"/><br>', plot_dir),
  sprintf('<h2>PCA</h2><img src="../%s/04_pca_treatment.png" width="600"/><br>', plot_dir),
  sprintf('<h2>UMAP</h2><img src="../%s/04_umap_treatment.png" width="600"/><br>', plot_dir),
  "</body></html>"
)

writeLines(html, report_file)
log_msg(sprintf("Generated HTML report at: %s", report_file))


cat("✅ Step 4 complete: Normalization, HVG selection, PCA & UMAP done and saved.\n")
