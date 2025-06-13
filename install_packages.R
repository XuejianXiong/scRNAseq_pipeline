# install_packages.R
# ------------------------
# Installs required R packages for single-cell RNA-seq analysis
# ------------------------

cran_packages <- c(
  "tidyverse",
  "Seurat",
  "patchwork",
  "Matrix",
  "data.table",
  "ggplot2",
  "readr",
  "dplyr",
  "cowplot",
  "sp",
  "shiny",
  "plotly",
  "reticulate",
  "stringi"
)

bioc_packages <- c(
  "SingleCellExperiment",
  "scater",
  "scran",
  "BiocManager",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "ComplexHeatmap"
)

# Install CRAN packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

cat("📦 Installing CRAN packages...\n")
invisible(lapply(cran_packages, install_if_missing))

# Install Bioconductor packages
cat("📦 Installing Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

cat("✅ All required R packages have been installed.\n")