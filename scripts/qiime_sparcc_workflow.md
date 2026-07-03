# QIIME2 + SCNIC Workflow for SparCC-based Network Analysis

This script outlines the steps used to generate SparCC-based correlation networks using the SCNIC plugin in QIIME 2.

### Requirements
- QIIME 2 (tested on v2021.8)
- SCNIC plugin installed (`qiime dev refresh-cache` required)
- Python 3.x for optional post-processing

### Steps
1. Generate BIOM and taxonomy files from phyloseq objects (R)
2. Convert BIOM files and import taxonomy into QIIME 2 (QIIME 2)
3. Filter low-abundance features and samples (QIIME 2)
4. Add pseudocounts and export the feature table (QIIME 2)
5. Calculate SparCC correlations and export correlation tables (SCNIC)
6. Calculate network topology and domain-specific metrics (Python)

See inline comments in this script for exact commands used in our pipeline.