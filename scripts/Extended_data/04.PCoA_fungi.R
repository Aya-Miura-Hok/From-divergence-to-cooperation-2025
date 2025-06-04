# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(microbiome)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(grid)

# Set ggplot theme
theme_set(theme_bw())

# Load data
sample_sheet <- read.csv("sample_data_env.csv", row.names = 1)
asv_sheet <- read.csv("1114otu_data.csv", row.names = 1)
tax_sheet <- read.csv("1114tax_data.csv", row.names = 1)

# Construct phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Subset bacterial taxa and sample types
ps<-subset_taxa(ps_all, Kingdom == "Fungi")
exp <- subset_samples(ps, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))


# Filter low-read samples
exp_filtered <- prune_samples(sample_sums(exp) > 966, exp)

# Collapse rare taxa and transform compositionally
pseq <- exp_filtered %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.1) %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Phylum") %>%
  microbiome::transform(transform = "compositional")

# Extract and merge OTU table and taxonomy
otu_table_raw <- t(as.data.frame(otu_table(pseq)))
taxonomy <- as.data.frame(tax_table(pseq)) %>% rownames_to_column("Sample")

otu_table_long <- otu_table_raw %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(taxonomy, by = "Sample")

# Extract environmental data
df_env <- as(sample_data(exp_filtered), "data.frame") %>%
  rownames_to_column("Sample") %>%
  select(Sample, Type, Site)

# Merge OTU and environment data
df_combined <- inner_join(otu_table_long, df_env, by = "Sample")

# Prepare OTU table for distance calculation
otu_table_numeric <- df_combined %>%
  select(-Sample, -Type, -Site) %>%
  mutate(across(everything(), as.numeric))

# Handle NAs and zero-sum rows
otu_table_numeric[is.na(otu_table_numeric)] <- 0
empty_rows <- rowSums(otu_table_numeric) == 0
otu_table_numeric <- otu_table_numeric[!empty_rows, ]
df_combined <- df_combined[!empty_rows, ]

# Compute Bray-Curtis distance
distance_matrix <- vegdist(otu_table_numeric, method = "bray")

# Run PCoA
pcoa_result <- cmdscale(distance_matrix, k = 2, eig = TRUE)
pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")

# Combine with metadata
sample_data_df <- df_combined %>% select(Sample, Type, Site)
pcoa_scores <- cbind(sample_data_df, pcoa_scores)

# Calculate vector fit (envfit)
envfit_result <- envfit(pcoa_result, otu_table_numeric, permutations = 999)
arrow_scores <- as.data.frame(scores(envfit_result, display = "vectors"))
colnames(arrow_scores) <- c("PCoA1", "PCoA2")
arrow_scores$Phylum <- rownames(arrow_scores)
arrow_scores$PCoA1 <- arrow_scores$PCoA1 * 0.2
arrow_scores$PCoA2 <- arrow_scores$PCoA2 * 0.2

# Plot PCoA with all vectors
pcoa_scores$Site <- factor(pcoa_scores$Site, levels = c("1","2","3","4"))
pcoa_scores$Type <- as.factor(pcoa_scores$Type)

p <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, color = Site, shape = Type)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Site, fill = Site), alpha = 0.2, geom = "polygon", 
               color = "black", linetype = "dashed", linewidth = 0.3) +
  theme_minimal() +
  labs(title = "PCoA with Phylum Contributions",
       subtitle = paste0("Explained Variance (PC1, PC2) = ", 
                         round((pcoa_result$eig[1] / sum(pcoa_result$eig)) * 100, 2), 
                         "%, ", 
                         round((pcoa_result$eig[2] / sum(pcoa_result$eig)) * 100, 2), 
                         "%")) +
  geom_segment(data = arrow_scores, aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.2, "cm")), 
               color = "red", linewidth = 1) +
  geom_text(data = arrow_scores, aes(x = PCoA1, y = PCoA2, label = Phylum), 
            inherit.aes = FALSE, color = "red", hjust = 1.2, size = 5)
print(p)

# Filter significant vectors
significant_phylum <- names(which(envfit_result$vectors$pvals < 0.05))
arrow_scores <- arrow_scores %>% filter(Phylum %in% significant_phylum)
label_offset <- 0.1
arrow_scores$label_x <- arrow_scores$PCoA1 + label_offset * sign(arrow_scores$PCoA1)
arrow_scores$label_y <- arrow_scores$PCoA2 + label_offset * sign(arrow_scores$PCoA2)

# Final PCoA plot with significant vectors only
p <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, color = Site, shape = Type)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Site, fill = Site), alpha = 0.2, geom = "polygon", 
               color = "black", linetype = "dashed", linewidth = 0.3) +
  theme_minimal() +
  labs(title = "PCoA with Significant Phylum Contributions",
       subtitle = paste0("Explained Variance (PC1, PC2) = ", 
                         round((pcoa_result$eig[1] / sum(pcoa_result$eig)) * 100, 2), 
                         "%, ", 
                         round((pcoa_result$eig[2] / sum(pcoa_result$eig)) * 100, 2), 
                         "%")) +
  geom_segment(data = arrow_scores, aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.3, "cm")), 
               color = "blue", linewidth = 1) +
  ggrepel::geom_text_repel(data = arrow_scores, aes(x = label_x, y = label_y, label = Phylum),
                           inherit.aes = FALSE, color = "blue", size = 5, 
                           segment.color = "blue", max.overlaps = 10) +
  theme(legend.position = "right",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14))
print(p)

# PERMANOVA and PERMDISP
permanova_result <- adonis2(distance_matrix ~ Site + Type, data = pcoa_scores, permutations = 2999)
dispersion_test <- betadisper(distance_matrix, pcoa_scores$Site)
dispersion_permutest <- permutest(dispersion_test, permutations = 2999)

print("PERMANOVA Result:")
print(permanova_result)

print("Homogeneity of Dispersions (PERMDISP) Result:")
print(dispersion_permutest)