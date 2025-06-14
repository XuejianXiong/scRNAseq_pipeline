# -----------------------------
# Step 5: Clustering and Marker Gene Analysis
# -----------------------------

# Clear environment, graphics, and console
rm(list = ls())            # Remove all variables from the workspace
graphics.off()             # Close all open graphics devices
cat("\014")                # Clear the console (works in RStudio)


library(Seurat)
library(ggplot2)
library(dplyr)

# -----------------------------
# Parameters and Paths
# -----------------------------

output_dir <- "results"
plot_dir <- "figures"
marker_dir <- file.path(output_dir, "markers")
input_file <- file.path(output_dir, "04_norm_dr_data.rds")
output_file <- file.path(output_dir, "05_clustered_data.rds")
report_file <- file.path(output_dir, "05_report.html")
log_file <- file.path(output_dir, "05_log.txt")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(marker_dir)) dir.create(marker_dir, recursive = TRUE)

params <- list(
  leiden_resolution = 0.5,
  n_neighbors = 15, 
  n_pcs = 40,
  top_n_genes = 10
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

log_msg("Step 5 started: Clustering and Marker Gene Analysis")
log_msg(sprintf(
  "Parameters - leiden_resolution: %.2f, n_neighbors: %d, n_pcs: %d, top_n_genes: %d",
  params$leiden_resolution, params$n_neighbors, params$n_pcs, params$top_n_genes
))

# -----------------------------
# Load Data
# -----------------------------

seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded normalized data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# Neighbors and Clustering
# -----------------------------

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:params$n_pcs, verbose = TRUE)
seurat_obj <- FindClusters(seurat_obj, resolution = params$leiden_resolution)

log_msg(sprintf("Identified clusters: %s", paste(levels(seurat_obj$seurat_clusters), collapse = ", ")))

# -----------------------------
# UMAP plots
# -----------------------------

seurat_obj <- RunUMAP(seurat_obj, dims = 1:params$n_pcs)

umap_cluster_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP: Clusters")
umap_cluster_file <- file.path(plot_dir, "05_umap_clusters.png")
ggsave(umap_cluster_file, umap_cluster_plot, width = 7, height = 5)

log_msg("Saved UMAP plot colored by clusters")

if ("treatment" %in% colnames(seurat_obj[[]])) {
  umap_treatment_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "treatment") + ggtitle("UMAP: Treatment")
  umap_treatment_file <- file.path(plot_dir, "05_umap_treatment.png")
  ggsave(umap_treatment_file, umap_treatment_plot, width = 7, height = 5)
  log_msg("Saved UMAP plot colored by treatment")
} else {
  log_msg("Metadata 'treatment' not found; skipping UMAP treatment plot")
}

# -----------------------------
# Marker Gene Identification
# -----------------------------

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

log_msg("Computed marker genes for each cluster")

# Save marker genes per cluster
clusters <- levels(seurat_obj$seurat_clusters)
summary_list <- list()

for (cl in clusters) {
  df <- markers %>%
    filter(cluster == cl) %>%
    select(cluster, gene, logfc = avg_log2FC, p_adj = p_val_adj)
  out_file <- file.path(marker_dir, sprintf("marker_genes_cluster_%s.csv", cl))
  write.csv(df, out_file, row.names = FALSE)
  log_msg(sprintf("Saved marker genes for cluster %s (%d genes)", cl, nrow(df)))
  summary_list[[cl]] <- head(df, params$top_n_genes)
}

# Save summary
summary_df <- bind_rows(summary_list)
write.csv(summary_df, "results/05_top_marker_genes_summary.csv", row.names = FALSE)
log_msg(sprintf("Saved top %d marker genes summary for all clusters", params$top_n_genes))

# -----------------------------
# HTML Report
# -----------------------------

html <- c(
  "<html><head><title>Step 5: Clustering & Marker Gene Analysis</title>",
  "<style>",
  "body { font-family: Arial, sans-serif; padding: 20px; }",
  "h1, h2 { color: #2c3e50; }",
  "img { border: 1px solid #ccc; padding: 10px; margin-bottom: 20px; width: 700px; }",
  "a { text-decoration: none; color: #2980b9; }",
  "</style></head><body>",
  "<h1>Step 5 Report: Leiden Clustering and Marker Genes</h1>",
  "<h2>UMAP Plots</h2>",
  sprintf('<img src="../%s/05_umap_clusters.png"><br>', plot_dir)
)

if ("treatment" %in% colnames(seurat_obj[[]])) {
  html <- c(html, sprintf('<img src="../%s/05_umap_treatment.png"><br>', plot_dir))
}

html <- c(html, "<h2>Marker Genes per Cluster</h2><ul>")
for (cl in clusters) {
  html <- c(html, sprintf('<li><a href="../%s/marker_genes_cluster_%s.csv">Marker genes for cluster %s</a></li>',
                          marker_dir, cl, cl))
}
html <- c(html,
          "</ul>",
          sprintf('<h2>Summary</h2><p><a href="../results/05_top_marker_genes_summary.csv">Top %d marker genes summary (CSV)</a></p>',
                  params$top_n_genes),
          "</body></html>"
)

writeLines(html, report_file)
log_msg(sprintf("Generated HTML report at %s", report_file))

# -----------------------------
# Save the data
# -----------------------------

saveRDS(seurat_obj, output_file)
log_msg(sprintf("Saved clustered Seurat object to %s", output_file))


cat("✅ Step 5 complete: Clustering and marker gene analysis done and saved.\n")
