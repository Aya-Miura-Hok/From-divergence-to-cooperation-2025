# QIIME2 + DADA2 Processing Workflow (Bacteria, Fungi, Eukaryote)

This document summarizes the preprocessing pipeline for converting raw FASTQ files into ASV tables using QIIME2 and DADA2. The procedure was applied to three datasets: **bacteria**, **fungi**, and **eukaryotes**.

---

## 1. Prepare Manifest Files

Create a manifest CSV file for each dataset (example: `manifest_bac.csv`, `manifest_fun.csv`, `manifest_euk.csv`) using absolute paths:

```csv
sample-id,absolute-filepath,direction
Sample01,/path/to/Sample01_R1.fastq.gz,forward
Sample01,/path/to/Sample01_R2.fastq.gz,reverse
...
```

## 2. Import FASTQ Files
```
# For bacteria
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_bac.csv \
  --output-path paired-end-demux_bac.qza \
  --input-format PairedEndFastqManifestPhred33

# For fungi
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_fun.csv \
  --output-path paired-end-demux_fun.qza \
  --input-format PairedEndFastqManifestPhred33

# For eukaryotes
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_euk.csv \
  --output-path paired-end-demux_euk.qza \
  --input-format PairedEndFastqManifestPhred33
```
