# Code for: "From divergence to cooperation: Microbial complementarity and symbiotic coexistence in stream biofilms"
This repository contains code used for statistical analyses and figure generation in the above manuscript, submitted to *Nature Ecology & Evolution*.

## Raw sequence data
The raw amplicon sequencing data (16S and 18S rRNA) have been deposited in the DNA Data Bank of Japan (DDBJ) under the BioProject accession number: **PRJDB16188**.

You can access the raw FASTQ files at:
[https://ddbj.nig.ac.jp/resource/bioproject/PRJDB16188](https://ddbj.nig.ac.jp/resource/bioproject/PRJDB16188)

Sample metadata and mapping files used for analysis are available in the `data/` folder of this repository.

## Folder structure
- `scripts/`: R scripts used for data processing, modeling, and visualization
- `data/`: Sample input data
- `scripts/setup.R`: Installs all required R packages

## Requirements
- R version 4.3.0 or later
- Required R packages are listed and installed via scripts/setup.R

## How to run
1. Clone this repository
Download or clone the repository to your local environment:
<pre><code> git clone https://github.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025.git 
 cd From-divergence-to-cooperation-2025 </code></pre>

2. Install R and required packages
Ensure you have R (version ≥ 4.3.0) installed.
Then, run the following in R or RStudio to install all required packages (including CRAN and Bioconductor packages):
<pre><code> source("scripts/setup.R") </code></pre>

3. Run analysis scripts
You can execute each script in the scripts/ folder depending on your interest (e.g., dbRDA, diversity metrics, pathway analysis).
For example:
<pre><code>source("scripts/dbRDA_analysis.R") </code></pre>

## Output
All plots and tables will be displayed directly in the R session.
No files will be saved to disk by default.
If you wish to export figures, please modify the relevant scripts to include ggsave() or similar output functions.

## Additional Workflow: QIIME2-based SparCC network analysis
This repository includes a reproducible pipeline for network analysis based on SparCC correlation metrics using QIIME 2 and SCNIC.

To run the QIIME2 workflow, please refer to the instructions in:
- `scripts/qiime_sparcc_workflow.md` (or `.sh` if you prefer shell scripts)

This workflow requires:
- QIIME 2 (tested with version 2021.8)
- SCNIC plugin for QIIME 2
- Python 3.x (for optional network formatting and hub detection)

Note: These steps assume you have already preprocessed the ASV table and taxonomy files using QIIME2-compatible formats. See the script for details.

## Preprocessing (FASTQ to ASV CSV)
Raw FASTQ files were processed using QIIME2 and DADA2 to generate ASV and taxonomy tables.  
All commands used for QIIME2 preprocessing are described in [`scripts/qiime_dada2_processing.md`](scripts/qiime_dada2_processing.md).
