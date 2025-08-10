# Required packages for the Upregulated Genes Pipeline
# Run this script to install all dependencies

# CRAN packages
cran_packages <- c(
  "ggplot2",
  "dplyr", 
  "enrichR"
)

# Bioconductor packages
bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

# Install CRAN packages
install.packages(cran_packages)

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(bioc_packages)

cat("✅ All required packages installed successfully!\n")
