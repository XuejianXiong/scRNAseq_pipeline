# -----------------------------
# Step 7: GSEA + Visualization + Report
# -----------------------------

# Clear environment, graphics, and console
rm(list = ls())            # Remove all variables from the workspace
graphics.off()             # Close all open graphics devices
cat("\014")                # Clear the console (works in RStudio)

library(Seurat)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(pheatmap)
library(knitr)
library(readr)
library(tibble)

source("R_scripts/pipeline_utils.R")

# -----------------------------
# Parameters and Paths
# -----------------------------
output_dir <- "results"
plot_dir <- "figures"
gsea_dir <- file.path(output_dir, "GSEA")
input_file <- file.path(output_dir, "06_de_data.rds")
html_report <- file.path(output_dir, "07_report.html")
log_file <- file.path(output_dir, "07_log.txt")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)

# Parameters for running GSEA with fgsea()
min_Size = 15
max_Size = 500
n_perm = 100


# -----------------------------
# Functions
# -----------------------------
# Function to NES barplot for top pathways
plot_barplot_nes <- function(df, title, save_path) {
  
  # Select top 10 pathways ranked by adjusted p-value (most significant)
  # Then sort them by NES for ordered plotting
  top <- df %>% 
    arrange(padj) %>%      # Sort by adjusted p-value (padj)
    head(10) %>%           # Take top 10 significant pathways
    arrange(NES)           # Sort these top 10 by NES for barplot order
  
  # Create barplot
  p <- ggplot(top, aes(
    NES,                                # X-axis: Normalized Enrichment Score
    reorder(pathway, NES),              # Y-axis: Pathways ordered by NES
    fill = NES                          # Fill color mapped to NES
  )) +
    geom_col() +                          # Draw bars (columns)
    scale_fill_gradient2(                 # Use diverging gradient: blue (low NES), white (mid), red (high NES)
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0                        # NES=0 is neutral enrichment
    ) +
    labs(
      title = title,                      # Plot title
      x = "Normalized Enrichment Score (NES)",  # X-axis label
      y = "Pathway"                       # Y-axis label
    ) +
    theme_minimal()                       # Clean minimal theme for clarity
  
  # Save plot to specified path
  ggsave(save_path, p, width = 8, height = 6, dpi = 150)
}


# Function to generate dotplot for top pathways
plot_dotplot_top_pathways <- function(df, title, save_path) {
  
  # Select top 10 pathways ranked by adjusted p-value (most significant)
  top <- df %>% arrange(padj) %>% head(10)
  
  # Scale gene set sizes relative to the largest size (for dot sizes)
  # Sizes are scaled to max 100 for easier interpretation
  sizes <- (top$size / max(top$size)) * 100
  
  # Create dotplot
  p <- ggplot(top, aes(
    NES,                                 # X-axis: Normalized enrichment score
    reorder(pathway, NES),               # Y-axis: Pathways ordered by NES
    size = sizes,                        # Dot size reflects gene set size
    color = -log10(padj)                 # Dot color reflects significance (-log10 of padj)
  )) +
    geom_point() +                         # Add the dot layer
    scale_color_viridis_c() +              # Use viridis color scale for good perceptual uniformity
    labs(
      title = title,                       # Plot title
      x = "NES",                           # X-axis label
      y = "Pathway",                       # Y-axis label
      size = "Gene Set Size",              # Legend for dot size
      color = "-log10(padj)"               # Legend for color scale
    ) +
    theme_minimal()                        # Use minimal theme for clean appearance
  
  # Save plot to specified path with given dimensions and resolution
  ggsave(save_path, p, width = 8, height = 6, dpi = 150)
}


# Function to generate heatmap of leading edge genes for the top pathway
plot_leading_edge_heatmap <- function(seurat_obj, lead_genes, group_by, group, save_path) {
  # Identify leading edge genes that are present in the Seurat object
  genes_in_data <- intersect(lead_genes, rownames(seurat_obj))
  
  # If no genes found, log message and exit the function early
  if (length(genes_in_data) == 0) {
    log_msg(sprintf("No leading edge genes found for %s=%s", group_by, group))
    return()
  }
  
  # Check which of the leading edge genes are already scaled
  scaled_genes <- rownames(GetAssayData(seurat_obj, slot = "scale.data"))
  missing_genes <- setdiff(genes_in_data, scaled_genes)
  
  # If any leading edge genes are missing from the scaled data, scale them
  if (length(missing_genes) > 0) {
    log_msg(sprintf("Scaling missing genes for %s=%s: %s", group_by, group, paste(missing_genes, collapse = ", ")))
    
    # Scale only the missing genes; returns a new Seurat object with scaled data for these genes
    scale_new <- Seurat::ScaleData(seurat_obj, features = missing_genes, verbose = FALSE)
    
    # Extract the existing scaled data and the new scaled data
    scale_data_full <- GetAssayData(seurat_obj, slot = "scale.data")
    scale_data_new <- GetAssayData(scale_new, slot = "scale.data")
    
    # Combine existing and new scaled data into a single matrix
    scale_data <- rbind(scale_data_full, scale_data_new)
    
    # Update the Seurat object with the combined scaled data matrix
    seurat_obj@assays[[DefaultAssay(seurat_obj)]]@scale.data <- scale_data
  }
  
  # Retrieve the (updated) scaled data matrix
  scale_data <- GetAssayData(seurat_obj, slot = "scale.data")
  
  # Subset the scaled data to just the leading edge genes
  expr <- scale_data[genes_in_data, , drop = FALSE]
  
  # Retrieve group metadata for ordering the columns (cells) in the heatmap
  meta <- seurat_obj@meta.data[[group_by]]
  
  # Order expression data columns by group metadata
  df_expr <- t(expr)[order(meta), ]
  
  # Generate and save the heatmap to file
  pheatmap(
    t(df_expr),                    # Transpose so genes are rows, cells are columns
    cluster_rows = FALSE,          # Do not cluster rows (genes)
    cluster_cols = FALSE,          # Do not cluster columns (cells)
    main = sprintf("Leading Edge Heatmap: %s=%s", group_by, group),  # Heatmap title
    filename = save_path,          # Path to save the heatmap image
    width = 8,                     # Width in inches
    height = 6                     # Height in inches
  )
}


log_msg("Step 7 started: GSEA and Visualization")

# -----------------------------
# Load data
# -----------------------------
seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded DE data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# GSEA (Gene Set Enrichment Analysis) on Seurat object
# -----------------------------
# Define grouping variables to perform GSEA across
groupings <- c("treatment", "seurat_clusters")

# Retrieve gene sets from MSigDB
# for both KEGG_LEGACY and KEGG_MEDICUS collections
gene_sets_legacy <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets_medicus <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")

# Combine both sets and structure as a named list of gene symbols by pathway
gene_sets <- bind_rows(gene_sets_legacy, gene_sets_medicus) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Loop through each grouping (e.g., treatment, cluster)
for (group_by in groupings) {
  
  # Skip if the grouping is not found in Seurat metadata
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    log_msg(sprintf("Skipping GSEA for '%s': metadata not found", group_by))
    next
  }
  
  # Identify unique groups within this grouping
  groups <- unique(seurat_obj@meta.data[[group_by]])
  
  # Perform GSEA for each group vs all others
  for (grp in groups) {
    log_msg(sprintf("Running GSEA for %s=%s", group_by, grp))
    
    # Differential expression: compare grp to all other cells
    # with no logFC threshold to retain all genes for ranking
    de <- FindMarkers(seurat_obj, ident.1 = grp, group.by = group_by, logfc.threshold = 0)
    
    # Create ranked list: rank genes by average log2 fold-change (descending)
    ranks <- de %>% 
      rownames_to_column("gene") %>% 
      arrange(desc(avg_log2FC)) %>% 
      select(gene, avg_log2FC) %>% 
      deframe() # Convert to named vector: names = gene, values = avg_log2FC
    
    # Run GSEA with fgsea, using specified min/max pathway sizes, and 100 permutations
    res <- fgsea(
      pathways = gene_sets, 
      stats = ranks, 
      minSize = min_Size, 
      maxSize = max_Size, 
      nperm = n_perm
      )

    # Sort results by adjusted p-value (padj)
    res <- res %>% arrange(padj)
    
    # Create output directory for results
    outdir <- file.path(gsea_dir, sprintf("%s_%s", group_by, grp))
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    # Save GSEA results as CSV
    csv_path <- file.path(outdir, "gsea_results.csv")
    write_csv(res, csv_path)
    log_msg(sprintf("Saved GSEA results to %s", csv_path))
    
    if (nrow(res) > 0) {
      # Extract leading edge genes from top enriched pathway (lowest padj)
      top_pathway <- res$leadingEdge[[1]]
      
      # Define output paths for plots
      bar_path <- file.path(outdir, sprintf("barplot_%s_%s.png", group_by, grp))
      dot_path <- file.path(outdir, sprintf("dotplot_%s_%s.png", group_by, grp))
      heatmap_path <- file.path(outdir, sprintf("heatmap_%s_%s.png", group_by, grp))
      
      # Generate NES barplot for top pathways
      plot_barplot_nes(res, sprintf("Top NES Barplot: %s=%s", group_by, grp), bar_path)
      log_msg(sprintf("Saved NES barplot for %s=%s at %s", group_by, grp, bar_path))
      
      # Generate dotplot for top pathways
      plot_dotplot_top_pathways(res, sprintf("Top Pathways Dotplot: %s=%s", group_by, grp), dot_path)
      log_msg(sprintf("Saved dotplot for %s=%s at %s", group_by, grp, dot_path))
      
      # Generate heatmap of leading edge genes for the top pathway
      plot_leading_edge_heatmap(seurat_obj, top_pathway, group_by, grp, heatmap_path)
      log_msg(sprintf("Saved leading edge heatmap for %s=%s at %s", group_by, grp, heatmap_path))
    } else {
      # No significant pathways found
      log_msg(sprintf("No significant enrichment for %s=%s", group_by, grp))
    }
  }
}

# -----------------------------
# HTML Report
# -----------------------------
# Initialize the HTML content as a character vector with header and report title
html <- c(
  # HTML header with title and main heading
  "<html><head><title>GSEA Report</title></head><body>",
  "<h1>GSEA and Visualization Report</h1>",
  # Add current timestamp into the report
  sprintf("<p><b>Date:</b> %s</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

# Loop through each grouping type (e.g., treatment, cluster)
for (group_by in groupings) {
  
  # Skip if the grouping column doesn't exist in metadata
  if (!(group_by %in% colnames(seurat_obj@meta.data))) next
  
  # Add a section heading for this grouping type
  html <- c(html, sprintf("<h2>GSEA Results by %s</h2>", tools::toTitleCase(group_by)))
  
  # Loop through each group (e.g., each treatment or cluster)
  groups <- unique(seurat_obj@meta.data[[group_by]])
  for (grp in groups) {
    
    # Define output directory and CSV path for this group
    outdir <- file.path("GSEA", sprintf("%s_%s", group_by, grp))
    csv_name <- file.path(outdir, "gsea_results.csv")
    
    # Load the top 10 GSEA results for preview
    df_preview <- read_csv(
      file.path(gsea_dir, sprintf("%s_%s", group_by, grp), "gsea_results.csv"),
      show_col_types = FALSE
    ) %>% head(10)
    
    # Convert the preview data frame to HTML table format using knitr::kable
    table_html <- knitr::kable(df_preview, format = "html", table.attr = "class='table table-striped'")
    
    # Append the section for this group: title, table, CSV link, and images
    html <- c(html, sprintf(
      "<h3>%s: %s</h3>
       <p>Top 10 GSEA results:</p>
       %s
       <p>Full CSV: <a href='%s'>%s</a></p>
       <img src='%s/barplot_%s_%s.png' width='600'><br>
       <img src='%s/dotplot_%s_%s.png' width='600'><br>
       <img src='%s/heatmap_%s_%s.png' width='600'><br><br>",
      tools::toTitleCase(group_by), grp,      # Group title
      table_html,                             # Top 10 results table
      csv_name, csv_name,                     # Link to full CSV
      outdir, group_by, grp,                  # Barplot image path
      outdir, group_by, grp,                  # Dotplot image path
      outdir, group_by, grp                   # Heatmap image path
    ))
  }
}

# Finalize the HTML by closing the tags
html <- c(html, "</body></html>")

# Write the HTML content to file
writeLines(html, html_report)
log_msg(sprintf("Generated HTML report at %s", html_report))


cat("✅ Step 7 complete: GSEA analysis and report generated.\n")