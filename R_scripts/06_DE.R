# -----------------------------
# Step 6: Differential Expression Analysis + Report
# -----------------------------

# Clear environment, graphics, and console
rm(list = ls())            # Remove all variables from the workspace
graphics.off()             # Close all open graphics devices
cat("\014")                # Clear the console (works in RStudio)

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------
# Parameters and Paths
# -----------------------------
output_dir <- "results"
plot_dir <- "figures"
csv_dir <- file.path(output_dir, "DE_csv")
input_file <- file.path(output_dir, "05_clustered_data.rds")
output_file <- file.path(output_dir, "06_de_data.rds")
html_report <- file.path(output_dir, "06_report.html")
log_file <- file.path(output_dir, "06_log.txt")

file_name1 <- "06_DE"
file_name2 <- "06_volcano"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

# Parameters for flagging genes as significant
cutoff_pval <- 0.05
cutoff_df <- 1


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


# Function to generate volcano plot from differential expression results
plot_volcano <- function(df, cutoff_pval, cutoff_df, title, save_path) {
  
  # Preprocess the input data frame:
  df <- df %>%
    mutate(
      # Replace NA adjusted p-values with 1 (non-significant)
      p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
      
      # Compute -log10 adjusted p-values for plotting
      # If p_val_adj is 0, replace with a small value to avoid -Inf
      log10_padj = -log10(ifelse(p_val_adj == 0, 1e-300, p_val_adj)),
      
      # Flag genes as significant if p_adj < 0.05 and abs(log2FC) > 1
      significant = (p_val_adj < cutoff_pval) & (abs(avg_log2FC) > cutoff_df)
    )
  
  # Create the volcano plot using ggplot2
  p <- ggplot(df, aes(x = avg_log2FC, y = log10_padj, color = significant)) +
    # Plot points, slightly transparent for better overlap visualization
    geom_point(alpha = 0.7) +
    
    # Color points: red = significant, grey = not significant
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    
    # Add horizontal line at significance threshold (-log10(0.05))
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    
    # Add vertical lines at log2 fold-change cutoffs (-1, 1)
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    
    # Add plot labels and theme
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value"
    ) +
    theme_minimal()
  
  # Save the plot to the specified file path
  ggsave(save_path, p, width = 8, height = 6, dpi = 150)
}


log_msg("Step 6 started: Differential Expression Analysis")

# -----------------------------
# Load data
# -----------------------------
# Load the pre-clustered Seurat object from the provided RDS file
seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded clustered data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# Differential Expression
# -----------------------------
# Define metadata columns to group by for DE analysis
groupings <- c("treatment", "seurat_clusters")

# Loop over each grouping variable (i.e. treatment, cluster)
for (group_by in groupings) {
  
  # Skip if the grouping column does not exist in Seurat metadata
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    log_msg(sprintf("Skipping DE for '%s': metadata not found", group_by))
    next
  }
  
  log_msg(sprintf("Starting DE analysis grouped by '%s'", group_by))
  
  # Run differential expression for all groups in the specified metadata column
  markers <- FindAllMarkers(seurat_obj, group.by = group_by, only.pos = FALSE)
  
  # Get the unique group labels (e.g., cluster IDs or treatments)
  unique_groups <- unique(markers$cluster)
  
  # Loop over each group to save DE results and generate plots
  for (grp in unique_groups) {
    
    # Filter markers for the current group
    df <- markers %>% filter(cluster == grp)
    
    # Save DE results as CSV
    csv_path <- file.path(csv_dir, sprintf(paste(file_name1, "_%s_%s.csv"), group_by, grp))
    write_csv(df, csv_path)
    log_msg(sprintf("Saved DE CSV for %s=%s at %s", group_by, grp, csv_path))
    
    # Generate and save volcano plot
    fig_path <- file.path(plot_dir, sprintf(paste(file_name2, "_%s_%s.png"), group_by, grp))
    fig_title <- sprintf("Volcano Plot: %s=%s", group_by, grp)
    
    # Call the volcano plot function to generate and save the plot
    plot_volcano(df, cutoff_pval, cutoff_df, fig_title, fig_path)
    log_msg(sprintf("Saved volcano plot for %s=%s at %s", group_by, grp, fig_path))
  }
}

# -----------------------------
# HTML Report
# -----------------------------
# Initialize the HTML content as a character vector with header and report title
html <- c(
  # HTML header with title and main heading
  "<html><head><title>Differential Expression Report</title></head><body>",
  "<h1>Differential Expression Analysis Report</h1>",
  # Add current timestamp into the report
  sprintf("<p><b>Date:</b> %s</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

# Loop over each grouping used for DE analysis (i.e., treatment, cluster)
for (group_by in groupings) {
  
  # Skip this grouping if it does not exist in Seurat metadata
  if (!(group_by %in% colnames(seurat_obj@meta.data))) next
  
  # Add section header for the grouping in the HTML report
  html <- c(html, sprintf("<h2>DE Results by %s</h2>", tools::toTitleCase(group_by)))
  
  # Loop over each unique group within this grouping
  for (grp in unique(seurat_obj@meta.data[[group_by]])) {
    
    # Build relative paths to CSV and volcano plot image
    csv_name <- sprintf(paste("DE_csv/", file_name1, "_%s_%s.csv"), group_by, grp)
    fig_name <- sprintf(paste("../figures/", file_name2, "_%s_%s.png"), group_by, grp)
    
    # Read the top 10 DE genes for preview in the report
    df_preview <- read_csv(file.path(csv_dir, sprintf(paste(file_name1, "_%s_%s.csv"), group_by, grp)), show_col_types = FALSE) %>%
      head(10)
    
    # Convert the preview dataframe to an HTML table with styling
    table_html <- knitr::kable(df_preview, format = "html", table.attr = "class='table table-striped'")
    
    # Append the DE results section for this group to the HTML content
    html <- c(html, sprintf(
      "<h3>%s: %s</h3>
       <p>Top 10 DE genes:</p>
       %s
       <p>Full CSV: <a href='%s'>%s</a></p>
       <img src='%s' alt='Volcano plot %s %s' width='600'><br><br>",
      tools::toTitleCase(group_by), grp,        # Section title (e.g., "Treatment: A")
      table_html,                               # Insert preview table
      csv_name, csv_name,                       # Link to full CSV file
      fig_name, group_by, grp                   # Insert volcano plot image
    ))
  }
}

# Close the HTML document
html <- c(html, "</body></html>")

# Write the HTML content to file
writeLines(html, html_report)
log_msg(sprintf("Generated HTML report at %s", html_report))

# -----------------------------
# Save Seurat object with DE results
# -----------------------------
saveRDS(seurat_obj, output_file)
log_msg(sprintf("Saved DE results to %s", output_file))


cat("✅ Step 6 complete: Differential expression analysis done and saved.\n")