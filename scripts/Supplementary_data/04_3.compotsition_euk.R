# ============================================================
# Input
# ============================================================
# This script requires a domain-specific phyloseq object
# generated in Script 00. asv_filter_euk.R.
#
# Available input file:
#   outputs/ps_final_euk.rds
#   outputs/ps_final_algae.rds
#   outputs/ps_final_protists.rds

# ============================================================
# Load libraries
# ============================================================
library(phyloseq)
library(dplyr)
library(ggplot2)

# Set default ggplot theme
theme_set(theme_bw())

# ============================================================
# Load phyloseq object
# ============================================================
dataset <- "euk"
# dataset <- "algae"
# dataset <- "protists"

ps_final <- readRDS(paste0("outputs/ps_final_", dataset, ".rds"))

# ============================================================
# Visualization (Phylum level_eukaryote)
# ============================================================
ps_f <- ps_final
tax <- tax_table(ps_f) %>% as.data.frame()
tax$Phylum <- ifelse(is.na(tax$Phylum) | tax$Phylum == "", "Unknown", tax$Phylum)
tax_table(ps_f) <- as.matrix(tax)

ps_phylum <- tax_glom(ps_f, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
df <- psmelt(ps_phylum_rel)

df$Phylum[df$Phylum %in% c("Unassigned", "unclassified", "unknown")] <- "Unknown"

top_phyla <- df %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

top_phyla <- unique(c(top_phyla, "Unknown"))

df <- df %>%
  mutate(Phylum = ifelse(Phylum %in% top_phyla, Phylum, "Other"))

df_summary <- df %>%
  group_by(Description, Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

df_summary <- df_summary %>%
  group_by(Description) %>%
  mutate(mean_abund = mean_abund / sum(mean_abund)) %>%
  ungroup()
  
df_summary$Phylum <- factor(df_summary$Phylum,
                            levels = c("Unknown", top_phyla[!top_phyla %in% c("Unknown")], "Other"))
							
euk_colors <- c(
  "Unknown"      = "#BFBFBF",
  "Other"        = "grey80",
  "Ciliophora"   = "#8491B4FF",
  "Chlorophyta"  = "#E69F00",
  "Diatomea"     = "#56B4E9",
  "Ochrophyta"   = "#0072B2",
  "Cercozoa"     = "#009E73",
  "Euglenozoa"   = "#CC79A7",
  "Apicomplexa"  = "#AD8BC9FF"
)

ggplot(df_summary, aes(x = Description, y = mean_abund, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = euk_colors, name = "Phylum") +
  labs(x = "Sample", y = "Relative abundance") +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.6),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
    axis.text.y = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.6, "cm")
  ) +
  theme(aspect.ratio = 0.8)

# ============================================================
# Visualization (Family level_algae group)
# ============================================================
ps_family <- tax_glom(ps_algae, taxrank = "Family")
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
ps_family_merged <- merge_samples(ps_family_rel, "Description")

ps_family_merged <- transform_sample_counts(ps_family_merged, function(x) x / sum(x))
df <- psmelt(ps_family_merged)
df$Family[df$Family %in% c("Unassigned", "unclassified", "unknown")] <- "Unknown"

top10_families <- df %>%
  group_by(Family) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abund)) %>%
  slice(1:10) %>%
  pull(Family)

df$Family[is.na(df$Family)] <- "Unknown"
df$Family[!df$Family %in% top10_families] <- "Other"
df$Family <- factor(df$Family, levels = c(top10_families, "Other"))

palette_algae <- c(
  "Bacillariophyceae" = "#DDAA66",
  "Ulvophyceae"       = "#6CBED5",
  "Sphaeropleales"    = "#66B2A3",
  "Chrysophyceae"     = "#5A6FAE",
  "Chlorophyceae"     = "#E59880",
  "Trebouxiophyceae"  = "#9B97B7",
  "Eustigmatales"     = "#8EC5BE",
  "Fragilariales"     = "#CC7B7B",
  "Chromulinales"     = "#9C8675",
  "Unknown"           = "#BFBFBF",
  "Other"             = "grey85"
)

base_levels <- c("Unknown", "Other", top10_families)
df$Family <- factor(df$Family, levels = unique(base_levels))

ggplot(df, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = palette_algae, name = "Family") +
  labs(x = "Sample", y = "Relative abundance") +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
    axis.text.y = element_text(size = 15)
  ) +
  theme(aspect.ratio = 0.8)

# ============================================================
# Visualization (Family level_protists group)
# ============================================================
ps_family <- tax_glom(ps_protists, taxrank = "Family")
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
ps_family_merged <- merge_samples(ps_family_rel, "Description")

ps_family_merged <- transform_sample_counts(ps_family_merged, function(x) x / sum(x))
df <- psmelt(ps_family_merged)
df$Family[df$Family %in% c("Unassigned", "unclassified", "unknown")] <- "Unknown"

top10_families <- df %>%
    group_by(Family) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice(1:10) %>%
    pull(Family)

df$Family <- as.character(df$Family)
df$Family[is.na(df$Family)] <- "Unknown"
df$Family[!df$Family %in% top10_families] <- "Other"
df$Family <- factor(df$Family, levels = c(top10_families, "Other"))

palette_protozoa <- c(
  "Oligohymenophorea" = "#DDAA66",
  "Eugregarinorida"   = "#6CBED5",
  "Heterotrichea"     = "#66B2A3",
  "Hypotrichia"       = "#5A6FAE",
  "Haptoria"          = "#E59880",
  "Phyllopharyngea"   = "#9B97B7",
  "Neobodonida"       = "#8EC5BE",
  "Colpodea"          = "#CC7B7B",
  "Rhizaspididae"     = "#9C8675",
  "Thecofilosea"      = "#9C8675",
  "Unknown"           = "#BFBFBF",
  "Other"             = "grey85"
)

base_levels <- c("Unknown", "Other", top10_families)
df$Family <- factor(df$Family, levels = unique(base_levels))

ggplot(df, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = palette_protozoa, name = "Family") +
  labs(x = "Sample", y = "Relative abundance") +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
    axis.text.y = element_text(size = 15)
  ) +
  theme(aspect.ratio = 0.8)