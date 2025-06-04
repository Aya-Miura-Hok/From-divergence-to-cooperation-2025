# QIIME2 + DADA2 Processing Workflow (Bacteria, Fungi, Eukaryote)

This document summarizes the preprocessing pipeline for converting raw FASTQ files into ASV tables using QIIME2 and DADA2. The procedure was applied to three datasets: **bacteria**, **fungi**, and **eukaryotes**.

---

## 1. Prepare Manifest Files

Create a manifest CSV file for each dataset (example: `manifest_bac.csv`, `manifest_fun.csv`, `manifest_euk.csv`) using absolute paths:

```csv
sample-id,absolute-filepath,direction
Sample01,/path/to/Sample01_R1.fastq.gz,forward
Sample01,/path/to/Sample01_R2.fastq.gz,reverse
... ```

## 1. Prepare Manifest Files

