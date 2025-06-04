# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Set default ggplot theme
theme_set(theme_bw())

# Load input data
data <- read_csv("e01.water_quality.csv")

# Filter rows with both δ13C and δ15N values available
data_filtered <- data %>%
  filter(!is.na(`δ13C`), !is.na(`δ15N`))

# Calculate mean and standard deviation of δ13C and δ15N by site and depth
summary_df <- data_filtered %>%
  group_by(Area, Water_Depth) %>%
  summarise(
    mean_d13C = mean(`δ13C`, na.rm = TRUE),
    sd_d13C = sd(`δ13C`, na.rm = TRUE),
    mean_d15N = mean(`δ15N`, na.rm = TRUE),
    sd_d15N = sd(`δ15N`, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    sd_d13C = replace_na(sd_d13C, 0),
    sd_d15N = replace_na(sd_d15N, 0),
    Site = as.factor(Area),
    Depth = factor(Water_Depth,
                   levels = c("surface water", "hyporheic 5cm water", "hyporheic 20cm water"))
  )

# Generate plot with error bars
ggplot(summary_df, aes(x = mean_d13C, y = mean_d15N)) +
  geom_point(aes(shape = Site, fill = Depth), size = 4, color = "black") +
  geom_errorbar(aes(ymin = mean_d15N - sd_d15N, ymax = mean_d15N + sd_d15N), width = 0.1) +
  geom_errorbarh(aes(xmin = mean_d13C - sd_d13C, xmax = mean_d13C + sd_d13C), height = 0.1) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +  # Distinguish sites by shape (fillable shapes)
  scale_fill_manual(
    values = c("black", "gray60", "white"),
    guide = guide_legend(override.aes = list(shape = 21))  # Ensure legend for fill uses fillable shape
  ) +
  labs(
    x = expression(delta^13*C~"(‰)"),
    y = expression(delta^15*N~"(‰)"),
    shape = "Site",
    fill = "Water Depth"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )