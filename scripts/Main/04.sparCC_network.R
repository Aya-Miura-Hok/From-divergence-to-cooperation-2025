# Reproducible workflow for SparCC-based network analysis of 16S data using QIIME2

# Load required libraries
library(phyloseq)
library(biomformat)
library(magrittr)

# --- [1] Import metadata, ASV table, and taxonomy ---
sample_sheet <- read.csv("sample_sheet_all.csv", row.names = 1)
asv_sheet <- read.csv("asv_sheet_all.csv", row.names = 1)
tax_sheet <- read.csv("tax_table_all.csv", row.names = 1)
asv_sheet <- t(asv_sheet)

ps <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# --- [2] Filter unassigned genera and export BIOM and taxonomy ---
taxonomy <- as.data.frame(tax_table(ps))
filtered_taxa <- rownames(taxonomy[!is.na(taxonomy$Genus), ])
ps_filtered <- prune_taxa(filtered_taxa, ps)
otu_table(ps_filtered) <- t(otu_table(ps_filtered))

otu_table(ps_filtered) %>%
  as("matrix") %>%
  make_biom() %>%
  write_biom("otu_table_genus.biom")

write.table(tax_table(ps_filtered), file = "taxonomy_genus.tsv", sep = "\t", quote = FALSE, col.names = NA)

metadata <- as(sample_data(ps), "data.frame")
metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]
write.table(metadata, file = "metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# --- [3] Import and format data in QIIME2 ---
# Convert BIOM to HDF5 format:
# biom convert -i otu_table_genus.biom -o otu_table_genus_v210.biom --table-type="OTU table" --to-hdf5

# Import feature table
# qiime tools import --input-path otu_table_genus_v210.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path feature-table.qza

# Format taxonomy file and import
# (Format Genus-level taxonomy into a "Feature ID" + "Taxon" structure with semicolon-separated prefixes)

# Export as "taxonomy_fixed.tsv", then:
# qiime tools import --input-path taxonomy_fixed.tsv --type 'FeatureData[Taxonomy]' --output-path taxonomy.qza

# If UTF-8 encoding fails:
# python3 -c 'with open("taxonomy_fixed.tsv", "r", encoding="ISO-8859-1") as infile, open("taxonomy_fixed_utf8.tsv", "w", encoding="utf-8") as outfile: [outfile.write(line) for line in infile]'

# --- [4] Feature/sample filtering ---
# Remove features with <50 reads:
# qiime feature-table filter-features --i-table feature-table.qza --p-min-frequency 50 --o-filtered-table filtered-table.qza

# Keep samples with >1000 reads:
# qiime feature-table filter-samples --i-table filtered-table.qza --p-min-frequency 1000 --o-filtered-table filtered-samples-table.qza

# Filter out rare ASVs (mean relative abundance < 0.01%)
# qiime feature-table filter-features --i-table filtered-samples-table.qza --p-min-frequency 10 --o-filtered-table filtered-abundance-table.qza

# --- [5] Add pseudocounts for compositional data ---
# qiime composition add-pseudocount --i-table filtered-abundance-table.qza --o-composition-table table-with-pseudocount.qza

# Convert to FeatureTable[Frequency] format (required for SCNIC)
# Export and re-import via:
# qiime tools export --input-path table-with-pseudocount.qza --output-path exported_data
# qiime tools import --type 'FeatureTable[Frequency]' --input-path exported_data/feature-table.biom --input-format BIOMV210Format --output-path table-frequency.qza

# Optional: log(x+1) transform for input table
# qiime tools export --input-path table-with-pseudocount.qza --output-path exported-frequency-table
# python log_transform.py
# biom convert -i log_transformed-table.tsv -o log_transformed-table.biom --table-type "OTU table" --to-hdf5
# qiime tools import --type 'FeatureTable[Frequency]' --input-path log_transformed-table.biom --output-path log-transformed-table.qza

# --- [6] SparCC network construction ---
# Calculate correlations using SparCC
# qiime SCNIC calculate-correlations --i-table table-frequency.qza --p-method sparcc --o-correlation-table correlations.qza

# Create modules based on positive correlations > 0.3
# qiime SCNIC make-modules-on-correlations \
#   --i-correlation-table correlations.qza \
#   --i-feature-table table-frequency.qza \
#   --p-min-r 0.3 \
#   --o-collapsed-table collapsed-table.qza \
#   --o-correlation-network correlation-network.qza \
#   --o-module-membership module-membership.qza

# --- [7] Export and visualize network ---
# Export correlation network to Cytoscape-compatible GML format:
# qiime tools export --input-path correlation-network.qza --output-path exported-network

# For signed networks (positive + negative edges):
# Use Python script (e.g., test3.py) to convert pairwise correlations to GML

# Identify hub/keystone species
# python keystones.py