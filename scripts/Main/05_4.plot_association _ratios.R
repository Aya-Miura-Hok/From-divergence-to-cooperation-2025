# ============================================================
# Input
# ============================================================
# Required input files in data/network_analysis/data/:
#   metrics_sites_signed.csv
#   domain_metrics_by_site_signed.csv

# ============================================================
# Load libraries
# ============================================================
library(tidyverse)
library(grid)

# ==============================================================================
# Plot Network Positive/Negative Association Ratios (overall and site-specific)
# ==============================================================================
site_csv <- "data/network_analysis/data/metrics_sites_signed.csv"
df_site  <- read.csv(site_csv, check.names = FALSE)

df_site_long <- df_site %>%
    rename(Site = 1) %>%
    select(Site, n_edges, pos_ratio, neg_ratio) %>%
    pivot_longer(cols = c(pos_ratio, neg_ratio),
                 names_to = "sign", values_to = "value") %>%
    mutate(
        sign = recode(sign, pos_ratio = "Positive", neg_ratio = "Negative"),
        Site = factor(Site, levels = c("All","Site1","Site2","Site3","Site4"))
    )

plot_donut_site <- function(df, group_cols = "Site") {
    df_sum <- df %>%
        group_by(across(all_of(group_cols)), sign) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
        group_by(across(all_of(group_cols))) %>%
        mutate(percent = value / sum(value) * 100)
    
    df_sum$sign <- factor(df_sum$sign, levels = c("Positive", "Negative"))
    
    ggplot(df_sum, aes(x = "", y = percent, fill = sign)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar(theta = "y", start = 0) +
        facet_wrap(as.formula(paste("~", group_cols)), nrow = 1) +
        scale_fill_manual(values = c(
            "Positive" = "#E59880",
            "Negative" = "#5A6FAE"
        )) +
        geom_text(
            aes(label = paste0(round(percent), "%")),
            position = position_stack(vjust = 0.5),
            size = 4,
            color = "white",
            family = "Helvetica",
            fontface = "bold"
        ) +
        theme_void(base_family = "Helvetica") +
        theme(
            strip.text = element_text(size = 13, face = "bold"),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
        )
}

p_site <- plot_donut_site(df_site_long, group_cols = "Site")
p_site

# ===================================================================
# Plot Network Positive/Negative Association Ratios (domain-specific)
# ===================================================================
dom_csv <- "data/network_analysis/data/domain_metrics_by_site_signed.csv"
df_dom  <- read.csv(dom_csv, check.names = FALSE)

names(df_dom)[1:2] <- c("Site", "Category")

df_dom_long <- df_dom %>%
  select(Site, Category, n_edges, pos_ratio, neg_ratio) %>%
  filter(!is.na(Category), Category != "") %>%
  mutate(
    Site = factor(Site, levels = c("All","Site1","Site2","Site3","Site4")),
    Category = factor(Category,
                      levels = c("Bac_Bac","Fun_Fun","Bac_Fun",
                                 "Bac_Euk","Fun_Euk","Euk_Euk"))
  ) %>%
  pivot_longer(cols = c(pos_ratio, neg_ratio),
               names_to = "sign", values_to = "value") %>%
  mutate(sign = recode(sign, pos_ratio = "Positive", neg_ratio = "Negative"))

plot_donut_domain <- function(df, group_cols = c("Site","Category")) {
  df_sum <- df %>%
    group_by(across(all_of(group_cols)), sign) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(percent = value / sum(value) * 100)

  df_sum$sign <- factor(df_sum$sign, levels = c("Positive", "Negative"))

  ggplot(df_sum, aes(x = "", y = percent, fill = sign)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y", start = 0) +
    facet_grid(Site ~ Category) +
    scale_fill_manual(values = c(
      "Positive" = "#E59880",
      "Negative" = "#5A6FAE"
    )) +
    geom_text(
      aes(label = paste0(round(percent), "%")),
      position = position_stack(vjust = 0.5),
      size = 3.2,
      color = "white",
      family = "Helvetica",
      fontface = "bold"
    ) +
    theme_void(base_family = "Helvetica") +
    theme(
      strip.text.x = element_text(size = 11, face = "bold"),
      strip.text.y = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
    )
}

p_dom <- plot_donut_domain(df_dom_long, group_cols = c("Site","Category")) +
    labs(y = "Site") +
    theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),
        strip.switch.pad.grid = unit(0.2, "cm")
    ) +
    facet_grid(Site ~ Category, switch = "y")
p_dom