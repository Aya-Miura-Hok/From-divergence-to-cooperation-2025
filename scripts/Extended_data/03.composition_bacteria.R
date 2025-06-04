# Load libraries
library(microViz)
library(microbiome)
library(phyloseq)
library(RColorBrewer)
theme_set(theme_bw())

# Load input data (metadata, OTU table, taxonomy)
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/metadata_bac.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/otu_table_exp.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/taxonomy_table.csv", row.names = 1)

# Create phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Subset bacterial taxa and relevant sample types
ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp <- subset_samples(ps, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))

# ---------- Phylum-level composition plot ---------- #
# Collapse rare taxa and aggregate by Description
pseq_phylum <- exp %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.45) %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Phylum") %>%
  merge_samples("Description") %>%
  microbiome::transform(transform = "compositional")

# Define 12-color palette
phylum_colors <- brewer.pal(12, "Paired")

# Plot Phylum-level composition
plot_phylum <- plot_composition(pseq_phylum, otu.sort = "abundance", transform = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(title = "Phylum"))

print(plot_phylum)

# ---------- Genus-level composition plot ---------- #
# Collapse rare taxa and aggregate by Description
pseq_genus <- exp %>%
  aggregate_rare(level = "Genus", detection = 5, prevalence = 0.55) %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Genus") %>%
  merge_samples("Description") %>%
  microbiome::transform(transform = "compositional")

# Define custom color palette for genus-level plot
genus_colors <- c(
  "#CAE0AB", "#7BAFDE", "#4EB265", "#F6C141", "#AA6F9E", "#DC050C",
  "#437DBF", "#CAACCB", "#F7F056", "#90C987", "#1965B0", "#D9CCE3",
  "#72190E", "#F1932D"
)

# Plot Genus-level composition
plot_genus <- plot_composition(pseq_genus, otu.sort = "abundance", transform = "identity") +
  scale_fill_manual(values = genus_colors) +
  guides(fill = guide_legend(title = "Genus"))

print(plot_genus)