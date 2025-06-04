# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyverse)
library(microbiome)
library(ggrepel)
library(RColorBrewer)

# Set default ggplot theme
theme_set(theme_bw())

# --- [1] Load metadata, OTU table, and taxonomy ---
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/metadata_bac.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/otu_table_exp.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/taxonomy_table.csv", row.names = 1)

ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# --- [2] Preprocessing: filtering and normalization ---
ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp <- subset_samples(ps, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))
exp_filtered <- prune_samples(sample_sums(exp) > 999, exp)

# --- [3] Normalize phylum-level abundance after filtering ---
pseq <- exp_filtered %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.45) %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Phylum") %>%
  microbiome::transform(transform = "compositional")

# --- [4] Construct feature matrix and merge environmental metadata ---
df_otu <- as(otu_table(pseq), "matrix") %>% t() %>%
  as.data.frame() %>% rownames_to_column("Sample")

df_env <- as(sample_data(pseq), "data.frame") %>%
  rownames_to_column("Sample") %>%
  select(Sample, Type, Site)

df_combined <- inner_join(df_otu, df_env, by = "Sample")

# --- [5] Load and normalize pathway abundance data ---
df_pathabun <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/pathabunｰfeature-table.biom.csv", row.names = 1, check.names = FALSE)
df_pathabun_norm <- df_pathabun / rowSums(df_pathabun)

# Select biologically relevant pathways (predefined)
selected_pathways <- c(
  "PWY-1042", "PWY-5100", "PWY-7456", "BENZOATE-COA-CAT-PWY",
  "PHENYLACETATE-CAT-PWY", "3-HYDROXYPHENYLACETATE-DEGRADATION-PWY", "GLYCOLYSIS",
  "PENTOSE-P-PWY", "ENTNER-DOUDOROFF-PWY", "ANAEROFRUCAT-PWY", "FASYN-INITIAL-PWY",
  "PWY0-1319", "PROPIONATE-PWY", "BUTANEDIOLFERMENTATION-PWY", "GLUTAMATE-DEG-PWY",
  "VALDEG-PWY", "ARGDEG-PWY", "TRNA-CHARGING-PWY", "PTS-PWY"
)

# Subset only pathways present in data
valid_pathways <- intersect(selected_pathways, colnames(df_pathabun_norm))

df_pathabun_filtered <- df_pathabun_norm %>%
  select(all_of(valid_pathways)) %>%
  rownames_to_column("Sample")

# Merge normalized pathways and phylum composition
df_combined_pwy <- inner_join(df_pathabun_filtered, df_combined, by = "Sample") %>%
  filter(complete.cases(.))

# --- [6] Format input matrices for dbRDA ---
df_X <- df_combined_pwy %>%
  select(Sample, Actinobacteriota:Verrucomicrobiota) %>%
  column_to_rownames("Sample")

df_Y <- df_combined_pwy %>%
  select(Sample, all_of(valid_pathways)) %>%
  column_to_rownames("Sample") %>%
  decostand("hellinger")  # Hellinger transformation

# --- [7] Run dbRDA and permutation tests ---
bray_dist <- vegdist(df_Y, method = "bray")
dbrda_model <- capscale(bray_dist ~ ., data = df_X)

anova_terms <- anova.cca(dbrda_model, by = "terms", permutations = 999)
sig_phylum <- rownames(anova_terms)[anova_terms$`Pr(>F)` < 0.05] |> na.omit()

# Re-run model with significant explanatory variables only
df_X_sig <- df_X[, sig_phylum]
common_samples <- intersect(rownames(df_Y), rownames(df_X_sig))
df_Y <- df_Y[common_samples, ]; df_X_sig <- df_X_sig[common_samples, ]

dbrda_model_new <- capscale(vegdist(df_Y, method = "bray") ~ ., data = df_X_sig)

# Inspect model summary
summary(dbrda_model_new)

# Global significance
anova.cca(dbrda_model_new, permutations = 999)

# Axis-level significance
anova.cca(dbrda_model_new, by = "axis", permutations = 999)

# Variable-level significance
anova.cca(dbrda_model_new, by = "terms", permutations = 999)

# --- [8] Extract sample and species scores ---
sample_scores <- scores(dbrda_model_new, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(df_combined_pwy %>% select(Sample, Site, Type), by = "Sample")

bp_scores <- scores(dbrda_model_new, display = "bp") %>% as.data.frame()

# Evaluate significance of pathway vectors
ko_vectors <- envfit(dbrda_model_new, df_Y, permutations = 999)
sig_pw <- names(ko_vectors$vectors$pvals[ko_vectors$vectors$pvals < 0.05])

species_scores <- as.data.frame(ko_vectors$vectors$arrows) %>%
  rownames_to_column("Pathway") %>%
  filter(Pathway %in% sig_pw)

species_scores$CAP1 <- species_scores$CAP1 * ko_vectors$vectors$r[sig_pw]
species_scores$CAP2 <- species_scores$CAP2 * ko_vectors$vectors$r[sig_pw]

# --- [9] Visualize dbRDA results (publication-ready figure) ---
ggplot() +
  geom_point(data = sample_scores,
             aes(x = CAP1, y = CAP2, fill = Site, shape = Type),
             size = 3.5, alpha = 0.9, color = "black") +

  geom_segment(data = bp_scores,
               aes(x = 0, y = 0, xend = CAP1 * 1.2, yend = CAP2 * 1.2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +

  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = CAP1 * 1.5, yend = CAP2 * 1.5),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "#F8766D", linetype = "dashed") +

  geom_text_repel(data = bp_scores,
                  aes(x = CAP1 * 1.7, y = CAP2 * 1.7, label = rownames(bp_scores)),
                  color = "black", size = 5) +

  geom_text_repel(data = species_scores,
                  aes(x = CAP1 * 2, y = CAP2 * 2, label = Pathway),
                  color = "#F8766D", size = 5) +

  scale_fill_brewer(palette = "Set2") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +

  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5, color = "black")),
    shape = guide_legend(override.aes = list(fill = "gray70"))
  ) +

  labs(
    title = "dbRDA: Significant Phyla and Pathways (P < 0.05)",
    x = "CAP1 (45.9%)", y = "CAP2 (24.2%)",
    fill = "Site", shape = "Type"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )
