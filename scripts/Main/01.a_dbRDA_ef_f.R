# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(microbiome)
library(tidyverse)
library(ggrepel)

# Set default ggplot theme
theme_set(theme_bw())

# Load input data (metadata, OTU table, taxonomy)
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/fungi/sample_data_env.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/fungi/1114otu_data.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/fungi/1114tax_data.csv", row.names = 1)

# Construct phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Subset bacterial taxa and select relevant sample types
ps<-subset_taxa(ps_all, Kingdom == "Fungi")
exp <- subset_samples(ps, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))

# Remove samples with fewer than 967 total reads
exp_filtered <- prune_samples(sample_sums(exp) > 966, exp)

# Collapse rare taxa at phylum level and perform compositional transformation
pseq <- exp_filtered %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.1) %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Phylum") %>%
  microbiome::transform(transform = "compositional")

# Extract OTU table and taxonomy and merge them
otu_table_raw <- t(as.data.frame(otu_table(pseq)))
taxonomy <- as.data.frame(tax_table(pseq)) %>% rownames_to_column("Sample")

otu_table_long <- otu_table_raw %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  left_join(taxonomy, by = "Sample")

# Extract environmental variables
df_env <- as(sample_data(exp_filtered), "data.frame") %>%
  rownames_to_column(var = "Sample") %>%
  select(Sample, Type, Site, d15N, d13C, HIX, BIX, POC, DOC)

# Merge OTU and environmental data
df_combined <- inner_join(otu_table_long, df_env, by = "Sample")

# Prepare response (phylum abundance) and explanatory (environmental) matrices
df_phylum <- df_combined %>%
  select(Sample, Ascomycota: Unknown) %>%
  column_to_rownames("Sample") %>%
  filter(complete.cases(.)) %>%
  filter(rowSums(.) > 0)

df_envs <- df_combined %>%
  select(Sample, d15N, d13C, HIX, BIX, POC, DOC) %>%
  column_to_rownames("Sample") %>%
  .[rownames(df_phylum), , drop = FALSE] %>%
  mutate_if(is.character, as.numeric)

# Apply Hellinger transformation and compute Bray-Curtis distances
df_phylum_hel <- decostand(df_phylum, method = "hellinger")
bray_dist <- vegdist(df_phylum_hel, method = "bray")
df_envs_scaled <- scale(df_envs)

# Perform distance-based redundancy analysis (dbRDA)
dbrda_model <- capscale(bray_dist ~ ., data = as.data.frame(df_envs_scaled))
summary(dbrda_model)

# Evaluate model and terms via permutation tests
anova.cca(dbrda_model, permutations = 999)
anova.cca(dbrda_model, by = "axis", permutations = 999)
anova.cca(dbrda_model, by = "terms", permutations = 999)

# Extract site scores and merge with sample metadata
sample_scores <- scores(dbrda_model, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(df_env[, c("Sample", "Type", "Site")], by = "Sample")
sample_scores$Site <- factor(as.character(sample_scores$Site), levels = c("1", "2", "3", "4"))
sample_scores$Type <- factor(sample_scores$Type, 
                             levels = c("leaf", "sediment_20", "sediment_5", "tile_d", "tile_l"))

# Fit phylum-level vectors to the ordination
phylum_vectors <- envfit(dbrda_model, df_phylum_hel, permutations = 999)
species_scores <- as.data.frame(phylum_vectors$vectors$arrows) %>%
  mutate(
    phylum = rownames(.),
    length = phylum_vectors$vectors$r,
    CAP1 = CAP1 * length,
    CAP2 = CAP2 * length
  )

# Identify significantly associated phyla (p < 0.05)
phylum_pvals <- data.frame(
  phylum = names(phylum_vectors$vectors$pvals),
  p_value = phylum_vectors$vectors$pvals
)

significant_phyla <- phylum_pvals %>%
  filter(p_value < 0.05) %>%
  pull(phylum)

species_scores_filtered <- species_scores %>%
  filter(phylum %in% significant_phyla)

# Extract scores of significant environmental variables
bp_scores <- scores(dbrda_model, display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  filter(Variable %in% c("d13C", "HIX", "BIX", "DOC"))

# Define color and shape palettes
site_colors <- c("1" = "#66C2A5", "2" = "#FC8D62", "3" = "#8DA0CB", "4" = "#E78AC3")
shape_values <- c("leaf" = 16, "sediment_20" = 17, "sediment_5" = 15, "tile_d" = 18, "tile_l" = 19)

# Generate final dbRDA plot
p <- ggplot() +
  geom_point(data = sample_scores, aes(x = CAP1, y = CAP2, color = Site, shape = Type),
             size = 3.2, alpha = 0.9, stroke = 0.6) +
  geom_segment(data = bp_scores, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue", linewidth = 1) +
  geom_text_repel(data = bp_scores, aes(x = CAP1 * 1.4, y = CAP2 * 1.4, label = Variable),
                  color = "blue", size = 5, segment.color = "grey50") +
  geom_segment(data = species_scores_filtered,
               aes(x = 0, y = 0, xend = CAP1 * 2, yend = CAP2 * 2, linewidth = length),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#8B0000", linetype = "dashed") +
  geom_text_repel(data = species_scores_filtered,
                  aes(x = CAP1 * 2.8, y = CAP2 * 2.8, label = phylum),
                  color = "#8B0000", size = 5,
                  box.padding = 0.5, point.padding = 0.8, segment.color = "grey50") +
  scale_color_manual(name = "Site", values = site_colors) +
  scale_shape_manual(name = "Type", values = shape_values) +
  scale_linewidth_continuous(name = "Impact (R^2)", range = c(0.3, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "dbRDA: Significant Phylum vs Environmental Factor (P < 0.05)",
    x = "CAP1 (51.1%)",
    y = "CAP2 (33.9%)"
  ) +
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 16, face = "bold")
  )

print(p)