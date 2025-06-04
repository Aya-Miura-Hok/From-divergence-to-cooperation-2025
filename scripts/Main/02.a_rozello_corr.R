# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(microbiome)
library(tidyverse)
library(ggrepel)
library(ggpmisc)

# Set default ggplot theme
theme_set(theme_bw())

# Load input data (metadata, OTU table, taxonomy)
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/eukaryote/sample_sheet_exp.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/eukaryote/otu_table_exp.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/eukaryote/taxonomy_table_exp.csv", row.names = 1)

# Construct phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Filter out samples with fewer than 1000 reads
read_counts <- sample_sums(ps_all)
exp_filtered <- prune_samples(read_counts > 999, ps_all)

# Compositional transformation (relative abundance)
pseq2 <- microbiome::transform(exp_filtered, transform = "compositional")

# Aggregate to phylum level
pseq3 <- aggregate_taxa(pseq2, level = "Phylum")

# Merge samples by site description
sample_data(pseq3)$Description <- factor(
  sample_data(pseq3)$Description,
  levels = sort(unique(sample_data(pseq3)$Description)),
  labels = paste0("Site_", sort(unique(sample_data(pseq3)$Description)))
)
pseq_grouped <- merge_samples(pseq3, "Description")
pseq_grouped <- microbiome::transform(pseq_grouped, transform = "compositional")

# Subset to 'Unknown' phylum only
Unknown <- subset_taxa(pseq_grouped, Phylum == "Unknown")

# Extract abundance data of 'Unknown' phylum
df_phylum <- psmelt(Unknown) %>%
  group_by(Sample, Phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Site = substr(Sample, nchar(Sample), nchar(Sample)))

# Boxplot to inspect 'Unknown' abundance by site
boxplot(Abundance ~ Site, data = df_phylum, main = "Unknown Abundance by Site")

# Load pre-processed Rozellomycota abundance data
df_phylum_Rozello <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/eukaryote/df_phylum_Rozello.csv", stringsAsFactors = FALSE)

# Merge 'Unknown' and Rozellomycota data by sample
merged_df <- left_join(df_phylum_Rozello, df_phylum, by = "Sample")

# Compute Pearson correlation
cor_result <- cor.test(
  merged_df$Abundance.x,
  merged_df$Abundance.y,
  method = "pearson"
)
print(cor_result)

# Plot correlation with regression line and statistics
ggplot(merged_df, aes(x = Abundance.x, y = Abundance.y)) +
  geom_point(color = "darkgreen", size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",
    label.y = "top",
    size = 6
  ) +
  labs(
    title = "Rozellomycota vs Unknown Abundance Correlation",
    x = "Rozellomycota Abundance",
    y = "Unknown Abundance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )
