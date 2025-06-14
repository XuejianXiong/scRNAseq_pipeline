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

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

# -----------------------------
# Logging Utility
# -----------------------------
log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n")
  write(line, file = log_file, append = TRUE)
}

log_msg("Step 6 started: Differential Expression Analysis")

# -----------------------------
# Volcano plot helper
# -----------------------------
plot_volcano <- function(df, title, save_path) {
  df <- df %>%
    mutate(
      p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
      log10_padj = -log10(ifelse(p_val_adj == 0, 1e-300, p_val_adj)),
      significant = (p_val_adj < 0.05) & (abs(avg_log2FC) > 1)
    )
  
  p <- ggplot(df, aes(x = avg_log2FC, y = log10_padj, color = significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
    theme_minimal()
  
  ggsave(save_path, p, width = 8, height = 6, dpi = 150)
}

# -----------------------------
# Load data
# -----------------------------
seurat_obj <- readRDS(input_file)
log_msg(sprintf("Loaded clustered data: %d cells × %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# Differential Expression
# -----------------------------
groupings <- c("treatment", "seurat_clusters")

for (group_by in groupings) {
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    log_msg(sprintf("Skipping DE for '%s': metadata not found", group_by))
    next
  }
  
  log_msg(sprintf("Starting DE analysis grouped by '%s'", group_by))
  markers <- FindAllMarkers(seurat_obj, group.by = group_by, only.pos = FALSE)
  
  unique_groups <- unique(markers$cluster)
  
  for (grp in unique_groups) {
    df <- markers %>% filter(cluster == grp)
    
    # Save CSV
    csv_path <- file.path(csv_dir, sprintf("06_DE_%s_%s.csv", group_by, grp))
    write_csv(df, csv_path)
    log_msg(sprintf("Saved DE CSV for %s=%s at %s", group_by, grp, csv_path))
    
    # Volcano plot
    fig_path <- file.path(plot_dir, sprintf("06_volcano_%s_%s.png", group_by, grp))
    plot_volcano(df, sprintf("Volcano Plot: %s=%s", group_by, grp), fig_path)
    log_msg(sprintf("Saved volcano plot for %s=%s at %s", group_by, grp, fig_path))
  }
}

# -----------------------------
# HTML Report
# -----------------------------
html <- c(
  "<html><head><title>Differential Expression Report</title></head><body>",
  "<h1>Differential Expression Analysis Report</h1>",
  sprintf("<p><b>Date:</b> %s</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

for (group_by in groupings) {
  if (!(group_by %in% colnames(seurat_obj@meta.data))) next
  
  html <- c(html, sprintf("<h2>DE Results by %s</h2>", tools::toTitleCase(group_by)))
  
  for (grp in unique(seurat_obj@meta.data[[group_by]])) {
    csv_name <- sprintf("DE_csv/06_DE_%s_%s.csv", group_by, grp)
    fig_name <- sprintf("../figures/06_volcano_%s_%s.png", group_by, grp)
    
    df_preview <- read_csv(file.path(csv_dir, sprintf("06_DE_%s_%s.csv", group_by, grp)), show_col_types = FALSE) %>%
      head(10)
    table_html <- knitr::kable(df_preview, format = "html", table.attr = "class='table table-striped'")
    
    html <- c(html, sprintf(
      "<h3>%s: %s</h3>
       <p>Top 10 DE genes:</p>
       %s
       <p>Full CSV: <a href='%s'>%s</a></p>
       <img src='%s' alt='Volcano plot %s %s' width='600'><br><br>",
      tools::toTitleCase(group_by), grp, table_html, csv_name, csv_name, fig_name, group_by, grp
    ))
  }
}

html <- c(html, "</body></html>")
writeLines(html, html_report)
log_msg(sprintf("Generated HTML report at %s", html_report))

# -----------------------------
# Save Seurat object
# -----------------------------
saveRDS(seurat_obj, output_file)
log_msg(sprintf("Saved DE results to %s", output_file))


cat("✅ Step 6 complete: Differential expression analysis done and saved.\n")
