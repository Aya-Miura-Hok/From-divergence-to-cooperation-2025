# Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)
library(RColorBrewer)

# Set default ggplot theme
theme_set(theme_bw())

# Load metadata, ASV table, and taxonomy (with FUNGuild annotations)
sample_sheet <- read.csv("sample_data_env.csv", row.names = 1)
asv_sheet <- read.csv("1114otu_data.csv", row.names = 1)
tax_sheet <- read.csv("funguild_input_guilds.csv", row.names = 1)

# Construct phyloseq object
ps_funguild <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Remove samples with fewer than 967 reads
read_counts_fun <- sample_sums(ps_funguild)
exp_filtered_fun <- prune_samples(read_counts_fun > 966, ps_funguild)

# Convert phyloseq object to long-format dataframe for plotting
df_tm <- psmelt(exp_filtered_fun)

# Clean Trophic Mode annotations
df_tm <- df_tm %>%
  mutate(
    Trophic.Mode = str_trim(Trophic.Mode),
    Trophic.Mode = ifelse(is.na(Trophic.Mode) | Trophic.Mode %in% c("-", "", NA),
                          "Unassigned",
                          Trophic.Mode)
  ) %>%
  mutate(Trophic.Mode = fct_collapse(as.factor(Trophic.Mode),
                                     Unassigned = c("Unassigned", "Unassigned ", " Unassigned", "-", "", NA)))

# Calculate relative abundance of each Trophic Mode per Site
site_summary <- df_tm %>%
  group_by(Site, Trophic.Mode) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Site) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))

# Define color palette for Trophic Modes
set3_colors <- brewer.pal(9, "Set3")
trophic_levels <- unique(site_summary$Trophic.Mode)
color_map <- setNames(set3_colors[1:length(trophic_levels)], trophic_levels)

# Plot: Trophic Mode composition
ggplot(site_summary, aes(x = Site, y = Relative_Abundance, fill = Trophic.Mode)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_map, name = "Trophic Mode") +
  labs(
    title = "Trophic Mode Composition",
    x = "Site",
    y = "Relative Abundance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )

 