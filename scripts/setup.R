# R version: 4.3.0 or later
# This script installs all required packages for running the analysis.

# 1. Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Define package groups
cran_packages <- c(
  "vegan", "ggplot2", "tidyverse", "ggrepel", "ggpmisc", "forcats", "stringr",
  "RColorBrewer", "magrittr", "readr", "tidyr", "car", "FSA", "rstatix",
  "grid", "iNEXT", "dplyr", "patchwork", "rlang"
)

bioc_packages <- c("phyloseq", "microbiome", "microViz", "biomformat")

# 3. Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# 4. Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

message("✔ All required packages have been checked/installed.")