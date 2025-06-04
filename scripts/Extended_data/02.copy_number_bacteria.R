# Load libraries
library(tidyverse)
library(phyloseq)
library(car)       # Levene test
library(FSA)       # Dunn’s test
library(rstatix)   # Kruskal-Wallis + Dunn

theme_set(theme_bw())

# Load data
sample_sheet <- read.csv("metadata_bac.csv", row.names = 1)
asv_sheet <- read.csv("otu_table_exp.csv", row.names = 1)
tax_sheet <- read.csv("taxonomy_table.csv", row.names = 1)

# Create phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Subset to Bacteria and target samples
ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp <- subset_samples(ps, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))

# ---------- Statistical Tests Per Substrate ---------- #
substrates <- c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20")

for (substrate in substrates) {
  cat("\n==== Substrate:", substrate, "====\n")
  
  ps_sub <- subset_samples(exp, Type == substrate)
  df_sub <- as.data.frame(as.matrix(sample_data(ps_sub)))
  
  # Convert and clean qPCR values
  df_sub$qpcr <- as.numeric(as.character(df_sub$qpcr))
  df_sub <- df_sub[!is.na(df_sub$qpcr), ]
  df_sub$log_qpcr <- log10(df_sub$qpcr)
  
  # Shapiro-Wilk test
  shapiro_test <- df_sub %>%
    group_by(Site) %>%
    summarise(p_value = shapiro.test(log_qpcr)$p.value)
  
  print("Shapiro-Wilk normality test:")
  print(shapiro_test)
  
  if (all(shapiro_test$p_value > 0.05)) {
    cat("→ ANOVA selected\n")
    result <- aov(log_qpcr ~ Site, data = df_sub)
    print(summary(result))
    print(TukeyHSD(result))
  } else {
    cat("→ Kruskal-Wallis + Dunn's test selected\n")
    print(kruskal.test(log_qpcr ~ Site, data = df_sub))
    print(
      df_sub %>% dunn_test(log_qpcr ~ Site, p.adjust.method = "bonferroni")
    )
  }
}

# ---------- Visualization ---------- #
# Full metadata extraction
df <- as.data.frame(as.matrix(sample_data(exp)))
df$qpcr <- as.numeric(as.character(df$qpcr))
df <- df[!is.na(df$qpcr), ]
df$log_qpcr <- log10(df$qpcr)

# Plot
p <- ggplot(df, aes(x = Site, y = log_qpcr, fill = Site)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  facet_wrap(~ Type, nrow = 1, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  labs(
    y = "log10 (copy number)",
    x = "Site"
  )

print(p)