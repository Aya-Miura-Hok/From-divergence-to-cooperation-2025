# QIIME2 + SCNIC Workflow for SparCC-based Network Analysis

This script outlines the steps used to generate correlation networks using the SCNIC plugin in QIIME 2 from 16S rRNA data.

### Requirements
- QIIME 2 (tested on v2024.2)
- SCNIC plugin installed (`qiime dev refresh-cache` required)
- Python 3.x for optional post-processing

### Steps
1. Convert BIOM + taxonomy to QIIME 2 input formats
2. Filter low-abundance features and samples
3. Add pseudocounts and export for SCNIC
4. Calculate SparCC correlations
5. Generate network modules and export GML
6. Identify keystone species via Python script

See inline comments in this script for exact commands used in our pipeline.
