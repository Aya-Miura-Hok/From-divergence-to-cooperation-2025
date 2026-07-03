# R version: 4.3.0 or later
# This script installs all required packages for running the analysis.

# 1. Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Define package groups
cran_packages <- c(
  "vegan", "ggplot2", "tidyverse", "ranger", "emmeans",
  "tidyr", "car", "FSA", "rstatix", "ape", 
  "iNEXT", "dplyr", "patchwork", "ggpubr"
)

bioc_packages <- c("phyloseq", "biomformat")

# 3. Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# 4. Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE, dependencies = TRUE)
  }
}

message("✔ All required packages have been checked/installed.")