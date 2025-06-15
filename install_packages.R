# install_packages.R
# ------------------------
# Installs required R packages for single-cell RNA-seq analysis
# ------------------------

# List of CRAN packages
cran_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "SeuratObject",
  "ggplot2",
  "patchwork",
  "readr",
  "tibble"
)

# Bioconductor packages
bioc_packages <- c(
  "fgsea",
  "msigdbr",
  "pheatmap",
  "knitr"
)

# -----------------------------
# Functions
# -----------------------------
# Function to install CRAN packages if missing
install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("📦 Installing CRAN package: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(sprintf("✔ CRAN package already installed: %s\n", pkg))
  }
}

# Function to install Bioconductor packages if missing
install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("📦 Installing Bioconductor package: %s\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("✔ Bioconductor package already installed: %s\n", pkg))
  }
}


# Install BiocManager first if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("📦 Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Deduplicate lists in case of overlap
cran_packages <- unique(cran_packages)
bioc_packages <- unique(bioc_packages)

# Install CRAN packages
cat("\n📌 Installing CRAN packages...\n")
invisible(lapply(cran_packages, install_if_missing_cran))

# Install Bioconductor packages
cat("\n📌 Installing Bioconductor packages...\n")
invisible(lapply(bioc_packages, install_if_missing_bioc))


cat("\n✅ All required R packages have been installed.\n")
