# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Load data
data <- read_csv("e01.water_quality.csv")

# Clean and reformat depth labels
data <- data %>%
  mutate(
    Water_Depth = trimws(Water_Depth),
    Depth = recode(Water_Depth,
                   "surface water" = "surface",
                   "hyporheic 5cm water" = "depth 5 cm",
                   "hyporheic 20cm water" = "depth 20 cm"),
    Site = as.character(Area),
    Group = paste(Site, Depth, sep = " - ")
  )

# Convert to long format for HIX and BIX indices
long_df <- data %>%
  pivot_longer(cols = c(HIX, BIX), names_to = "Index", values_to = "Value") %>%
  drop_na(Value)

# Calculate mean and standard deviation by group and index
summary_df <- long_df %>%
  group_by(Group, Index) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    .groups = "drop"
  )

# Define desired order of Site and Depth for plotting
site_order <- c("1", "2", "3", "4")
depth_order <- c("surface", "depth 5 cm", "depth 20 cm")
group_levels <- as.vector(t(outer(site_order, depth_order, paste, sep = " - ")))

# Reorder factor levels for proper axis arrangement
summary_df$Group <- factor(summary_df$Group, levels = group_levels)

# Pivot wider for plotting
wide_df <- summary_df %>%
  pivot_wider(names_from = Index, values_from = c(Mean, SD)) %>%
  mutate(Group = factor(Group, levels = group_levels))

# Generate dual-axis plot with legend
ggplot(wide_df, aes(x = Group)) +
  geom_point(aes(y = Mean_HIX, color = "HIX"), shape = 21, fill = "black", size = 3) +
  geom_errorbar(aes(ymin = Mean_HIX - SD_HIX, ymax = Mean_HIX + SD_HIX, color = "HIX"), width = 0.2) +
  geom_point(aes(y = Mean_BIX * 6, color = "BIX"), shape = 22, fill = "gray60", size = 3) +
  geom_errorbar(aes(ymin = (Mean_BIX - SD_BIX) * 6, ymax = (Mean_BIX + SD_BIX) * 6, color = "BIX"), width = 0.2) +
  scale_color_manual(
    name = "",
    values = c("HIX" = "black", "BIX" = "gray60")
  ) +
  scale_y_continuous(
    name = "HIX",
    sec.axis = sec_axis(~./6, name = "BIX")
  ) +
  labs(x = "Site and Depth") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 13),
    axis.title.y.right = element_text(size = 13, color = "gray30"),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )