# Load required libraries
library(phyloseq)
library(iNEXT)

# Load data
sample_sheet <- read.csv("metadata_bac.csv", row.names = 1)
asv_sheet <- read.csv("otu_table_exp.csv", row.names = 1)
tax_sheet <- read.csv("taxonomy_table.csv", row.names = 1)

# Construct phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Subset for Bacteria and target sample types
ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp <- subset_samples(ps, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))
exp_filtered <- prune_samples(sample_sums(exp) > 999, exp)

# Create list of phyloseq objects per site
site_names <- unique(sample_data(exp_filtered)$Site)
site_list <- lapply(site_names, function(site) subset_samples(exp_filtered, Site == site))
names(site_list) <- site_names

# Convert to iNEXT input list by site
otu_list_by_site <- lapply(site_list, function(physeq_site) {
  otu_mat <- as(otu_table(physeq_site), "matrix")
  otu_mat <- otu_mat[rowSums(otu_mat) > 1, , drop = FALSE]
  split(t(otu_mat), rownames(otu_mat))
})

# Run iNEXT for each site
iNEXT_results <- lapply(seq_along(site_names), function(i) {
  cat("Running iNEXT for site:", site_names[i], "\n")
  iNEXT(otu_list_by_site[[i]], q = c(0, 1, 2), datatype = "abundance", endpoint = 1000)
})
names(iNEXT_results) <- site_names

# Plot curves
par(mfrow = c(2, 2))
for (i in seq_along(iNEXT_results)) {
  plot(iNEXT_results[[i]], type = 1, main = paste("Site:", site_names[i]))
}

# Collect all results
collect_results <- function(results_list, type) {
  all_results <- list()
  all_cols <- character()
  
  for (i in seq_along(site_names)) {
    res <- results_list[[i]]$iNextEst[[type]]
    if (!is.null(res)) {
      res$Site <- site_names[i]
      all_cols <- union(all_cols, colnames(res))
      all_results[[i]] <- res
    }
  }
  
  all_results <- lapply(all_results, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (m in missing) df[[m]] <- NA
    df[, all_cols]
  })
  
  do.call(rbind, all_results)
}

size_df <- collect_results(iNEXT_results, "size_based")
coverage_df <- collect_results(iNEXT_results, "coverage_based")

# Save to CSV
write.csv(size_df, "iNEXT_Sitewise_Results_size_based_site.csv", row.names = FALSE)
write.csv(coverage_df, "iNEXT_Sitewise_Results_coverage_based_site.csv", row.names = FALSE)

# Load required libraries
library(dplyr)         # Data manipulation
library(ggplot2)       # Plotting
library(FSA)           # For Dunn's test
library(ggsignif)      # Significance annotation on plots
library(tidyr)         # Data reshaping
library(purrr)         # For map function
library(gridExtra)     # To arrange multiple plots

# Load iNEXT results
coverage_based_results <- read.csv("iNEXT_Sitewise_Results_coverage_based_site.csv")

# Filter rows with coverage ~0.95
coverage_95_results <- dplyr::filter(coverage_based_results, abs(SC - 0.95) < 0.01)

# Extract diversity indices at q = 0, 1, 2
diversity_at_95_coverage <- coverage_95_results %>%
  dplyr::filter(Order.q %in% c(0, 1, 2)) %>%
  dplyr::select(Site, Order.q, qD, qD.LCL, qD.UCL)

# Kruskal-Wallis and Dunn’s test with annotated boxplots
kruskal_results <- list()
dunn_results <- list()
plots <- list()

for (q in c(0, 1, 2)) {
  cat("\n-----------------------------\n")
  cat("Order.q = ", q, "\n")
  cat("-----------------------------\n")

  data_q <- dplyr::filter(diversity_at_95_coverage, Order.q == q)
  data_q$Site <- as.factor(data_q$Site)

  # Kruskal-Wallis test
  kruskal_test <- kruskal.test(qD ~ Site, data = data_q)
  kruskal_results[[paste0("Order.q_", q)]] <- kruskal_test
  print(kruskal_test)

  # Dunn's test
  dunn_test <- dunnTest(qD ~ Site, data = data_q, method = "bonferroni")
  dunn_results[[paste0("Order.q_", q)]] <- dunn_test

  # Handle column names
  colnames_lower <- tolower(colnames(dunn_test$res))
  comp_col <- grep("comparison", colnames_lower, value = TRUE)
  padj_col <- grep("p.adj", colnames_lower, value = TRUE)

  # Extract significant pairs
  significant_pairs <- dunn_test$res %>%
    dplyr::filter(!!as.name(padj_col) < 0.05) %>%
    dplyr::select(all_of(comp_col), all_of(padj_col))

  # Annotate if significant
  if (nrow(significant_pairs) > 0) {
    signif_annotations <- significant_pairs %>%
      tidyr::separate(1, into = c("Site1", "Site2"), sep = " - ") %>%
      mutate(
        Site1 = as.character(Site1),
        Site2 = as.character(Site2),
        y_position = max(data_q$qD, na.rm = TRUE) + (1:n()) * (0.1 * max(data_q$qD, na.rm = TRUE)),
        annotations = if_else(!!as.name(padj_col) < 0.01, "**", "*")
      )

    comparisons_list <- signif_annotations %>%
      dplyr::select(Site1, Site2) %>%
      split(1:nrow(signif_annotations)) %>%
      purrr::map(~ as.character(c(.x$Site1, .x$Site2)))

    plot <- ggplot(data_q, aes(x = Site, y = qD, fill = Site)) +
      geom_boxplot() +
      labs(title = paste0("Order.q = ", q), x = "Site", y = "qD (Diversity Index)") +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      geom_signif(
        comparisons = comparisons_list,
        annotations = signif_annotations$annotations,
        y_position = signif_annotations$y_position,
        tip_length = 0.005,
        textsize = 6
      )
  } else {
    plot <- ggplot(data_q, aes(x = Site, y = qD, fill = Site)) +
      geom_boxplot() +
      labs(title = paste0("Order.q = ", q), x = "Site", y = "qD (Diversity Index)") +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
  }

  plots[[paste0("Order.q_", q)]] <- plot
}

# Combine all plots
gridExtra::grid.arrange(
  plots$Order.q_0,
  plots$Order.q_1,
  plots$Order.q_2,
  nrow = 1
)
