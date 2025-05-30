# Code for: "From divergence to cooperation: Microbial complementarity and symbiotic coexistence in stream biofilms"
This repository contains code used for statistical analyses and figure generation in the above manuscript, submitted to *Nature Ecology & Evolution*.

## Raw sequence data
The raw amplicon sequencing data (16S and 18S rRNA) have been deposited in the DNA Data Bank of Japan (DDBJ) under the BioProject accession number: **PRJDB16188**.

You can access the raw FASTQ files at:
[https://ddbj.nig.ac.jp/resource/bioproject/PRJDB16188](https://ddbj.nig.ac.jp/resource/bioproject/PRJDB16188)

Sample metadata and mapping files used for analysis are available in the `data/` folder of this repository.

## Folder structure
- `scripts/`: R scripts used for data processing, modeling, and visualization
- `data/`: Sample input data (for demonstration; full dataset available upon request or via Zenodo)

## Requirements
- R version 4.3.0 or later
- Packages: vegan, picante, ggplot2, etc. (see `scripts/setup.R`)

## How to run
1. Clone this repository
2. Run `scripts/main_analysis.R`
3. Outputs will be saved in `output/`

## License
This code is licensed under the MIT License. See LICENSE file for details.
