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
library(tidyverse)
library(grid)

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
# dbRDA analysis
# ============================================================
# Aggregate ASVs at the Family level
ps_family <- tax_glom(ps_final, taxrank = "Family", NArm = TRUE)
ps_family <- prune_taxa(taxa_sums(ps_family) > 0, ps_family)

# Create a Hellinger-transformed Community Matrix
otu <- as(otu_table(ps_family), "matrix")
if (!taxa_are_rows(otu_table(ps_family))) otu <- t(otu)

rel <- sweep(otu, 2, colSums(otu), "/")
rel[!is.finite(rel)] <- 0

hell <- sqrt(rel)

# Match Metadata to Community Data
meta <- data.frame(sample_data(ps_final))
samps <- colnames(hell)
meta  <- meta[rownames(meta) %in% samps, , drop = FALSE]
meta  <- meta[match(samps, rownames(meta)), , drop = FALSE]
stopifnot(identical(rownames(meta), samps))

hell_use <- t(hell)
meta_use <- meta

# Prepare and Standardize Environmental Variables
env_candidates <- c("d13C", "DOC", "POC", "HIX", "BIX", "TDN", "TDP")
avail <- intersect(env_candidates, colnames(meta_use))

if (length(avail) == 0) {
  stop("No environmental variables were found in the metadata. Please check column names.")
}

env_use <- meta_use[, avail, drop = FALSE] |>
    lapply(function(x) as.numeric(as.character(x))) |>
    as.data.frame()

# Remove samples with missing environmental data
keep_rows <- complete.cases(env_use)
hell_use  <- hell_use[keep_rows, , drop = FALSE]
env_use   <- env_use[keep_rows, , drop = FALSE]
meta_use  <- meta_use[keep_rows, , drop = FALSE]

# Standardize environmental variables
env_std <- as.data.frame(scale(env_use))

# ============================================================
# Fit the Full Partial dbRDA Model
# ============================================================
rhs <- paste(colnames(env_std), collapse = " + ")
form_txt <- sprintf("hell_use ~ %s + Condition(replicate)", rhs)

df_all   <- cbind(env_std, replicate = meta_use$replicate)

cap_full <- capscale(as.formula(form_txt), data = df_all, distance = "bray")

# ============================================================
# Evaluate the Full Model
# ============================================================
# Overall test
print(anova(cap_full, permutations = 999, strata = meta_use$replicate))

# Term-wise test
print(anova(cap_full, permutations = 999, by = "terms", strata = meta_use$replicate))

# Adjusted R-squared
print(RsquareAdj(cap_full))

# Variables with VIF > 10 may indicate strong multicollinearity
print(vif.cca(cap_full))

# ============================================================
# Forward Selection Based on Adjusted R-squared
# ============================================================
cap_null  <- capscale(hell_use ~ 1 + Condition(replicate), data = df_all, distance = "bray")

set.seed(123)
cap_sel <- ordiR2step(cap_null, scope = formula(cap_full),
                      permutations = 999, trace = 1)

cap_final <- cap_sel

# ============================================================
# Evaluate the Final dbRDA Model
# ============================================================
print(anova(cap_final, permutations = 999, strata = meta_use$replicate))
print(anova(cap_final, permutations = 999, by = "terms", strata = meta_use$replicate))
print(RsquareAdj(cap_final))

# ============================================================
# Plot dbRDA Results
# ============================================================
scr_sites <- scores(cap_final, display = "sites", choices = 1:2) |> as.data.frame()
scr_bipl  <- scores(cap_final, display = "bp",    choices = 1:2) |> as.data.frame()

plot_df <- scr_sites |>
    rownames_to_column("SampleID") |>
    left_join(meta_use |> rownames_to_column("SampleID"), by = "SampleID") |>
    mutate(Site = as.factor(Site),
           Type = as.factor(Type))

arrow_mul <- 1
scr_bipl2 <- scr_bipl |> 
    mutate(ar_x = CAP1 * arrow_mul, ar_y = CAP2 * arrow_mul,
           var  = rownames(scr_bipl))

eig   <- cap_final$CCA$eig
prop1 <- eig[1] / sum(eig)
prop2 <- eig[2] / sum(eig)

p <- ggplot(plot_df, aes(CAP1, CAP2, color = Site, shape = Type)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Site), linetype = 2, level = 0.68, show.legend = FALSE) +
  geom_segment(
    data = scr_bipl2,
    aes(x = 0, y = 0, xend = ar_x, yend = ar_y),
    arrow = arrow(length = unit(0.18, "cm")),
    inherit.aes = FALSE, color = "black", linewidth = 0.4
  ) +
  geom_text(
    data = scr_bipl2,
    aes(x = ar_x, y = ar_y, label = var),
    size = 4, vjust = -0.4, family = "Arial", fontface = "plain",
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    "1" = "#D55E00",  # reddish brown
    "2" = "#009E73",  # green
    "3" = "#0072B2",  # blue
    "4" = "#CC79A7"   # purple-pink
  )) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    aspect.ratio = 1,
	axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = sprintf("dbRDA1 (%.1f%%)", 100 * prop1),
    y = sprintf("dbRDA2 (%.1f%%)", 100 * prop2)
  )

print(p)