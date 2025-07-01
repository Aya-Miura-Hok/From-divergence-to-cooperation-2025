# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rlang)

# Set global theme
theme_set(theme_bw())

# Load data
diversity <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/iNEXT_Sitewise_Results_coverage_based_site.csv")
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/From-divergence-to-cooperation-2025/refs/heads/main/data/bacteria/metadata_bac.csv", row.names = 1)
sample_sheet$SampleName <- rownames(sample_sheet)

# Extract coverage ~0.95 and q = 0, 1, 2
diversity_95 <- diversity %>%
  filter(abs(SC - 0.95) < 0.01 & Order.q %in% c(0, 1, 2)) %>%
  rename(SampleName = Assemblage)

# Merge with metadata
merged_df <- merge(diversity_95, sample_sheet, by = "SampleName")

# Convert environmental variables to numeric
for (var in c("HIX", "d13C", "DOC", "POC")) {
  merged_df[[var]] <- as.numeric(as.character(merged_df[[var]]))
}

# Generic plotting function for correlation between qD and environmental variables
plot_qD_vs_env <- function(df, q_val, var_name, x_label, cor_position = "topright") {
  df_q <- df %>% filter(Order.q == q_val)
  cor_result <- cor.test(df_q$qD, df_q[[var_name]], method = "spearman")
  
  x_pos <- if (cor_position == "topleft") min(df_q[[var_name]]) else max(df_q[[var_name]])
  y_pos <- max(df_q$qD)
  hjust_val <- if (cor_position == "topleft") 0 else 1
  
  ggplot(df_q, aes(x = !!sym(var_name), y = qD)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
    labs(
      title = paste("q =", q_val),
      x = x_label,
      y = "Diversity (qD)"
    ) +
    annotate("text",
             x = x_pos, y = y_pos,
             hjust = hjust_val, vjust = 1.2,
             label = paste0("r = ", round(cor_result$estimate, 2),
                            "\nP = ", signif(cor_result$p.value, 3)),
             size = 5
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
}

# Plot for each environmental variable
plot_env_group <- function(var_name, x_label, cor_position = "topright") {
  p0 <- plot_qD_vs_env(merged_df, 0, var_name, x_label, cor_position)
  p1 <- plot_qD_vs_env(merged_df, 1, var_name, x_label, cor_position)
  p2 <- plot_qD_vs_env(merged_df, 2, var_name, x_label, cor_position)
  p0 + p1 + p2 + plot_layout(ncol = 3)
}

# Generate plots
p1 <- plot_env_group("HIX", "HIX")
p2 <- plot_env_group("d13C", "δ13C")
p3 <- plot_env_group("DOC", "DOC", cor_position = "topleft")
p4 <- plot_env_group("POC", "POC")

print(p1)
print(p2)
print(p3)
print(p4)