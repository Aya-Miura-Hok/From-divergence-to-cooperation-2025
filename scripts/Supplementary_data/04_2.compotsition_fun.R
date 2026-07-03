# ============================================================
# Input
# ============================================================
# This script requires a domain-specific phyloseq object
# generated in Script 00. asv_filter_fun.R.
#
# Available input file:
#   outputs/ps_final_fun.rds

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
ps_final <- readRDS("outputs/ps_final_fun.rds")

# ============================================================
# Visualization (Phylum level)
# ============================================================
ps_f <- ps_final
tax <- tax_table(ps_f) %>% as.data.frame()
tax$Phylum <- ifelse(is.na(tax$Phylum) | tax$Phylum == "", "Unknown", tax$Phylum)
tax_table(ps_f) <- as.matrix(tax)

ps_phylum <- tax_glom(ps_f, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
df <- psmelt(ps_phylum_rel)

top_phyla <- df %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 5) %>%
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
                    levels = c("Unknown", "Ascomycota", "Basidiomycota",
                               "Chytridiomycota", "Rozellomycota", "Other"))

phylum_colors <- c(
  "Ascomycota" = "#8491B4FF",
  "Basidiomycota" = "#E69F00",
  "Chytridiomycota" = "#56B4E9",
  "Rozellomycota" = "#CC79A7",
  "Unknown" = "#BFBFBF",
  "Other" = "grey80"
)

ggplot(df_summary, aes(x = Description, y = mean_abund, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = phylum_colors) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
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
# Visualization (Family level)
# ============================================================
ps_family <- tax_glom(ps_f, taxrank = "Family")
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
ps_family_merged <- merge_samples(ps_family_rel, "Description")

ps_family_merged <- transform_sample_counts(ps_family_merged, function(x) x / sum(x))
df <- psmelt(ps_family_merged)

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

palette_pastel <- c(
  "Rozellomycota_fam_Incertae_sedis" = "#DDAA66",
  "Cladosporiaceae"                 = "#6CBED5",
  "Mycosphaerellaceae"              = "#66B2A3",
  "Hyaloscyphaceae"                 = "#5A6FAE",
  "Discinellaceae"                  = "#E59880",
  "Leotiaceae"                      = "#9B97B7",
  "Halosphaeriaceae"                = "#8EC5BE",
  "GS08_fam_Incertae_sedis"         = "#CC7B7B",
  "Saccotheciaceae"                 = "#9C8675",
  "Helotiaceae"                     = "#B9A998",
  "Other"                           = "grey85"
)

df$Family <- factor(df$Family, levels = c("Other", setdiff(levels(df$Family), "Other")))

ggplot(df, aes(x = Sample, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = palette_pastel, name = "Family") +
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
    ) + theme(aspect.ratio = 0.8)