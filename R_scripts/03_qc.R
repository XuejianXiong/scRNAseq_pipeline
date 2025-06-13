# -----------------------------
# Step 3: Quality Control
# -----------------------------

library(Seurat)
library(ggplot2)
library(dplyr)

# -----------------------------
# Parameters and Paths
# -----------------------------

input_file <- "results/02_merged_data.rds"
output_file <- "results/03_filtered_data.rds"
qc_plot_dir <- "figures"
report_file <- "results/03_report_R.html"
log_file <- "results/03_log.txt"

qc_params <- list(
  min_genes = 200,
  max_pct_mt = 10,
  min_cells_per_gene = 3
)

if (!dir.exists(qc_plot_dir)) dir.create(qc_plot_dir, recursive = TRUE)

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
log_msg(sprintf("Loaded merged data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# QC Metrics
# -----------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# -----------------------------
# QC Plot Helper
# -----------------------------

save_qc_plots <- function(seurat_obj, group_label, suffix) {
  vln_file <- file.path(qc_plot_dir, sprintf("qc_violin_by_%s_%s.png", group_label, suffix))
  scatter_mt_file <- file.path(qc_plot_dir, sprintf("qc_scatter_mt_by_%s_%s.png", group_label, suffix))
  scatter_gene_file <- file.path(qc_plot_dir, sprintf("qc_scatter_genes_by_%s_%s.png", group_label, suffix))
  
  g1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = group_label, pt.size = 0.1, ncol = 3)
  ggsave(vln_file, g1, width = 10, height = 4)
  
  g2 <- FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt")
  ggsave(scatter_mt_file, g2, width = 10, height = 5)
  
  g3 <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA")
  ggsave(scatter_gene_file, g3, width = 10, height = 5)
  
  log_msg(sprintf("Saved QC plots (%s) for %s", suffix, group_label))
}

# -----------------------------
# QC Before Filtering
# -----------------------------

group_cols <- c("treatment")

for (group in group_cols) {
  save_qc_plots(seurat_obj, group, "before")
}

# -----------------------------
# Filtering
# -----------------------------

pre_shape <- dim(seurat_obj)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= qc_params$min_genes & percent.mt <= qc_params$max_pct_mt
)

# Filter genes based on min_cells_per_gene
gene_counts <- Matrix::rowSums(GetAssayData(seurat_obj, assay = "RNA", layer = "counts") > 0)
genes_to_keep <- names(gene_counts[gene_counts >= qc_params$min_cells_per_gene])
seurat_obj <- subset(seurat_obj, features = genes_to_keep)

post_shape <- dim(seurat_obj)

log_msg(sprintf("Filtering complete — Before: %d cells × %d genes; After: %d cells × %d genes",
                pre_shape[2], pre_shape[1], post_shape[2], post_shape[1]))

log_msg(sprintf("Applied thresholds — min_genes: %d, min_cells_per_gene: %d, max_pct_mito: %d%%",
                qc_params$min_genes, qc_params$min_cells_per_gene, qc_params$max_pct_mt))

# -----------------------------
# QC After Filtering
# -----------------------------

for (group in group_cols) {
  save_qc_plots(seurat_obj, group, "after")
}

# -----------------------------
# Save Filtered Data
# -----------------------------

saveRDS(seurat_obj, output_file)
log_msg(sprintf("Filtered data saved to %s", output_file))

# -----------------------------
# Simple HTML Report
# -----------------------------

html <- c(
  "<html><head><title>QC Report</title></head><body>",
  "<h1>Single-Cell QC Report</h1>",
  sprintf("<p><b>Date:</b> %s</p>", Sys.time()),
  sprintf("<p><b>Initial shape:</b> %d cells × %d genes</p>", pre_shape[2], pre_shape[1]),
  sprintf("<p><b>Post-filtering shape:</b> %d cells × %d genes</p>", post_shape[2], post_shape[1]),
  sprintf("<p><b>QC thresholds:</b> min_genes = %d, min_cells_per_gene = %d, max_pct_mito = %d%%</p>",
          qc_params$min_genes, qc_params$min_cells_per_gene, qc_params$max_pct_mt),
  "<h2>QC Metrics Before Filtering</h2>"
)

for (group in group_cols) {
  html <- c(html,
            sprintf('<img src="../%s/qc_violin_by_%s_before.png" width="600"/><br>', qc_plot_dir, group),
            sprintf('<img src="../%s/qc_scatter_mt_by_%s_before.png" width="600"/><br>', qc_plot_dir, group),
            sprintf('<img src="../%s/qc_scatter_genes_by_%s_before.png" width="600"/><br>', qc_plot_dir, group))
}

html <- c(html, "<h2>QC Metrics After Filtering</h2>")

for (group in group_cols) {
  html <- c(html,
            sprintf('<img src="../%s/qc_violin_by_%s_after.png" width="600"/><br>', qc_plot_dir, group),
            sprintf('<img src="../%s/qc_scatter_mt_by_%s_after.png" width="600"/><br>', qc_plot_dir, group),
            sprintf('<img src="../%s/qc_scatter_genes_by_%s_after.png" width="600"/><br>', qc_plot_dir, group))
}

html <- c(html, "</body></html>")
writeLines(html, report_file)
log_msg(sprintf("HTML QC report generated at: %s", report_file))

cat("✅ Step 3 complete: Quality control applied, filtered data saved.\n")
