# ============================================================
# Input
# ============================================================
# This script requires a domain-specific phyloseq object
# generated in Script 00_1. asv_filter_bac.R, 00_2. asv_filter_fun.R, 00_3. asv_filter_euk.R.
#
# Available input files:
#   outputs/ps_final_bac.rds
#   outputs/ps_final_fun.rds
#   outputs/ps_final_algae.rds
#   outputs/ps_final_protists.rds
#
# Select the dataset to analyze by loading the corresponding RDS file.

# ============================================================
# Load libraries
# ============================================================
library(phyloseq)
library(vegan)
library(ape)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Set default ggplot theme
theme_set(theme_bw())

# ============================================================
# Load phyloseq object
# ============================================================
dataset <- "bac"
# dataset <- "fun"
# dataset <- "algae"
# dataset <- "protists"

ps_final <- readRDS(paste0("outputs/ps_final_", dataset, ".rds"))

# ============================================================
# PCoA analysis
# ============================================================
# Aggregate ASVs at the Family level
ps_family <- tax_glom(ps_final, taxrank = "Family")

# Extract OTU table and apply Hellinger transformation
otu <- as(otu_table(ps_family), "matrix")
if(!taxa_are_rows(otu_table(ps_family))) otu <- t(otu)

rel <- sweep(otu, 2, colSums(otu), "/"); rel[!is.finite(rel)] <- 0
hell <- sqrt(rel)

# Calculate Bray–Curtis dissimilarity and perform PCoA
bray <- vegdist(t(hell), method = "bray")
pcoa <- ape::pcoa(bray)

# Prepare metadata for visualization
meta <- data.frame(sample_data(ps_final))
scores <- data.frame(pcoa$vectors[,1:2], SampleID = rownames(pcoa$vectors))
colnames(scores)[1:2] <- c("PCoA1","PCoA2")
plot_df <- scores %>%
  left_join(meta %>% tibble::rownames_to_column("SampleID"), by = "SampleID")
plot_df$Site <- factor(plot_df$Site)

# Check homogeneity of multivariate dispersion before PERMANOVA
bd_site <- betadisper(bray, meta$Site)
anova(bd_site)

# ============================================================
# PERMANOVA
# ============================================================
# Perform PERMANOVA with replicate as a blocking factor
perm <- adonis2(bray ~ Site + Type, data = meta, 
        strata = meta$replicate, permutations = 999)
print(perm)

# Optional: test interaction if sample size allows
# perm2 <- adonis2(bray ~ Site * Type, data = meta,
#                  strata = meta$replicate, permutations = 999)
# print(perm2)

# Estimate marginal effect sizes (pseudo-R²)
perm_terms <- vegan::adonis2(
  bray ~ Site + Type,
  data = meta,
  strata = meta$replicate,
  permutations = 999,
  by = "margin"
)

perm_terms_df <- as.data.frame(perm_terms) |>
  tibble::rownames_to_column("term") |>
  dplyr::filter(term %in% c("Site","Type")) |>
  dplyr::select(term, Df, SumOfSqs, R2, F = `F`, p.value = `Pr(>F)`)

print(perm_terms_df)

# ============================================================
# PCoA plot
# ============================================================
gg <- ggplot(plot_df, aes(PCoA1, PCoA2, color = Site, shape = Type)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Site), linetype = 2, level = 0.68, show.legend = FALSE) +
  scale_color_manual(values = c(
    "1" = "#D55E00",
    "2" = "#009E73",
    "3" = "#0072B2",
    "4" = "#CC79A7"
  )) +
  theme_classic(base_size = 16, base_family = "Arial") +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", 100 * pcoa$values$Relative_eig[1]),
    y = sprintf("PCoA2 (%.1f%%)", 100 * pcoa$values$Relative_eig[2])
  ) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    aspect.ratio = 1,
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

print(gg)

if (!dir.exists("outputs")) dir.create("outputs")

write.csv(perm_terms_df,
          paste0("outputs/permanova_family_", dataset, ".csv"),
          row.names = FALSE)

ggsave(
    filename = "Fig3A.pdf",
    plot = gg,
    device = cairo_pdf,
	width = 180,
    height = 180,
    units = "mm",
    bg = "white"
)