# -----------------------------
# Step 3: Quality Control
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
input_file <- file.path(output_dir, "02_merged_data.rds")
output_file <- file.path(output_dir, "03_filtered_data.rds")
report_file <- file.path(output_dir, "03_report.html")
log_file <- file.path(output_dir, "03_log.txt")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

qc_params <- list(
  min_genes = 200,
  max_pct_mt = 10,
  min_cells_per_gene = 3
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

# Function to generate and save QC plots
save_qc_plots <- function(seurat_obj, group_label, suffix) {
  # Construct file paths for each plot 
  # (violin, scatter with percent.mt, scatter with nFeature_RNA)
  vln_file <- file.path(plot_dir, sprintf("qc_violin_by_%s_%s.png", group_label, suffix))
  scatter_mt_file <- file.path(plot_dir, sprintf("qc_scatter_mt_by_%s_%s.png", group_label, suffix))
  scatter_gene_file <- file.path(plot_dir, sprintf("qc_scatter_genes_by_%s_%s.png", group_label, suffix))
  
  # Generate violin plots of key QC metrics 
  # (number of features, counts, % mitochondrial content)
  g1 <- VlnPlot(
    seurat_obj, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = group_label,    # Group cells by the specified metadata (e.g., treatment, sample)
    pt.size = 0.1,             # Make individual cell points small to reduce overplotting
    ncol = 3                   # Arrange violin plots in 3 columns
  )
  # Save the violin plot
  ggsave(vln_file, g1, width = 10, height = 4)  
  
  # Generate scatter plot: total counts vs. % mitochondrial reads
  g2 <- FeatureScatter(
    seurat_obj, 
    "nCount_RNA", 
    "percent.mt"
  )
  # Save the scatter plot
  ggsave(scatter_mt_file, g2, width = 10, height = 5)  
  
  # Generate scatter plot: total counts vs. number of detected genes
  g3 <- FeatureScatter(
    seurat_obj, 
    "nCount_RNA", 
    "nFeature_RNA"
  )
  # Save the scatter plot
  ggsave(scatter_gene_file, g3, width = 10, height = 5) 
  
  log_msg(sprintf("Saved QC plots (%s) for %s", suffix, group_label))
}


log_msg("Step 3 started: Quality Control")

# -----------------------------
# Load Data
# -----------------------------
seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded merged data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# QC Metrics
# -----------------------------
# Calculate the percentage of reads mapping to mitochondrial genes (genes starting with "MT-")
# and add it as a metadata column named 'percent.mt'
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# List of metadata columns to group QC plots
group_cols <- c("treatment") 

# Generate and save QC plots for each group before filtering
for (group in group_cols) {
  save_qc_plots(seurat_obj, group, "before")
}

# -----------------------------
# Filtering
# -----------------------------
# Record the initial dimensions (genes × cells) before filtering
pre_shape <- dim(seurat_obj)

# Filter out low-quality cells:
# - Keep cells with at least min_genes expressed
# - Keep cells with mitochondrial percentage below the threshold
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= qc_params$min_genes & percent.mt <= qc_params$max_pct_mt
)

# Filter out lowly expressed genes:
# - Count how many cells each gene is detected in (non-zero counts)
gene_counts <- Matrix::rowSums(GetAssayData(seurat_obj, assay = "RNA", layer = "counts") > 0)

# - Retain genes expressed in at least min_cells_per_gene cells
genes_to_keep <- names(gene_counts[gene_counts >= qc_params$min_cells_per_gene])

# - Subset Seurat object to keep only those genes
seurat_obj <- subset(seurat_obj, features = genes_to_keep)

# Record dimensions after filtering
post_shape <- dim(seurat_obj)

# Log the filtering summary (cells and genes before vs. after)
log_msg(sprintf(
  "Filtering complete — Before: %d cells × %d genes; After: %d cells × %d genes",
  pre_shape[2], pre_shape[1], post_shape[2], post_shape[1]
))

# Log applied QC thresholds for traceability
log_msg(sprintf(
  "Applied thresholds — min_genes: %d, min_cells_per_gene: %d, max_pct_mito: %d%%",
  qc_params$min_genes, qc_params$min_cells_per_gene, qc_params$max_pct_mt
))

# -----------------------------
# QC After Filtering
# -----------------------------
# Generate and save QC plots for each group after filtering
for (group in group_cols) {
  save_qc_plots(seurat_obj, group, "after")
}

# Save the filtered Seurat object to file for downstream analysis
saveRDS(seurat_obj, output_file)
log_msg(sprintf("Filtered data saved to %s", output_file))


# -----------------------------
# HTML Report
# -----------------------------
# Initialize HTML content as a character vector (each element = one line)
html <- c(
  "<html><head><title>QC Report</title></head><body>",
  "<h1>Single-Cell QC Report</h1>",
  # Insert the current date and time
  sprintf("<p><b>Date:</b> %s</p>", Sys.time()),
  # Insert initial data shape before filtering (cells × genes)
  sprintf("<p><b>Initial shape:</b> %d cells × %d genes</p>", pre_shape[2], pre_shape[1]),
  # Insert data shape after filtering
  sprintf("<p><b>Post-filtering shape:</b> %d cells × %d genes</p>", post_shape[2], post_shape[1]),
  # Insert QC thresholds used for filtering
  sprintf("<p><b>QC thresholds:</b> min_genes = %d, min_cells_per_gene = %d, max_pct_mito = %d%%</p>",
          qc_params$min_genes, qc_params$min_cells_per_gene, qc_params$max_pct_mt),
  "<h2>QC Metrics Before Filtering</h2>"
)

# Add plots for each group (e.g. "treatment") before filtering
for (group in group_cols) {
  html <- c(html,
            # Violin plot of QC metrics
            sprintf('<img src="../%s/qc_violin_by_%s_before.png" width="600"/><br>', plot_dir, group),
            # Scatter plot: nCount_RNA vs percent.mt
            sprintf('<img src="../%s/qc_scatter_mt_by_%s_before.png" width="600"/><br>', plot_dir, group),
            # Scatter plot: nCount_RNA vs nFeature_RNA
            sprintf('<img src="../%s/qc_scatter_genes_by_%s_before.png" width="600"/><br>', plot_dir, group))
}

# Add section header for after-filtering plots
html <- c(html, "<h2>QC Metrics After Filtering</h2>")

# Add plots for each group after filtering
for (group in group_cols) {
  html <- c(html,
            # Violin plot of QC metrics
            sprintf('<img src="../%s/qc_violin_by_%s_after.png" width="600"/><br>', plot_dir, group),
            # Scatter plot: nCount_RNA vs percent.mt
            sprintf('<img src="../%s/qc_scatter_mt_by_%s_after.png" width="600"/><br>', plot_dir, group),
            # Scatter plot: nCount_RNA vs nFeature_RNA
            sprintf('<img src="../%s/qc_scatter_genes_by_%s_after.png" width="600"/><br>', plot_dir, group))
}

# Close the HTML tags
html <- c(html, "</body></html>")

# Write the HTML content to file
writeLines(html, report_file)
log_msg(sprintf("HTML QC report generated at: %s", report_file))


cat("✅ Step 3 complete: Quality control applied, filtered data saved.\n")
