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
library(ranger)

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
# Random Forest analysis
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

# Prepare Environmental Variables
env_vars <- c("d13C", "DOC", "POC", "HIX", "BIX", "TDN", "TDP")
avail <- intersect(env_vars, colnames(meta_use))

if (length(avail) == 0) {
  stop("No environmental variables were found in the metadata. Please check column names.")
}

env_use <- meta_use[, avail, drop = FALSE] |>
  lapply(function(x) as.numeric(as.character(x))) |>
  as.data.frame()

rownames(env_use) <- rownames(meta_use)

# Remove samples with missing environmental data
keep_rows <- complete.cases(env_use)
hell_use  <- hell_use[keep_rows, , drop = FALSE]
env_use   <- env_use[keep_rows, , drop = FALSE]
meta_use  <- meta_use[keep_rows, , drop = FALSE]

# Standardize environmental variables
env_std <- as.data.frame(scale(env_use))
stopifnot(identical(rownames(env_std), rownames(meta_use)))

# ============================================================
# Calculate PCoA Axes from Community Composition
# ============================================================
pc <- cmdscale(vegdist(hell_use, "bray"), k = 2, eig = TRUE)

pc_df <- tibble(
  SampleID = rownames(hell_use),
  PC1 = pc$points[,1],
  PC2 = pc$points[,2]
) %>%
  left_join(env_std |> rownames_to_column("SampleID"), by="SampleID", suffix=c("", "")) %>%
  left_join(meta_use |> rownames_to_column("SampleID"), by="SampleID", suffix=c("", ""))

# ============================================================
# Run Random Forest Models
# ============================================================
dat  <- pc_df |>
  dplyr::select(all_of(c("PC1","PC2", avail))) |>
  tidyr::drop_na()

set.seed(123)
rf_pc1 <- ranger::ranger(
  PC1 ~ ., data = dplyr::select(dat, PC1, all_of(avail)),
  num.trees = 2000, importance = "permutation", oob.error = TRUE
)
rf_pc2 <- ranger::ranger(
  PC2 ~ ., data = dplyr::select(dat, PC2, all_of(avail)),
  num.trees = 2000, importance = "permutation", oob.error = TRUE
)

# ============================================================
# Summarize Model Performance and Variable Importance
# ============================================================
oob_r2 <- function(mod, y) 1 - mod$prediction.error / stats::var(y)
cat(sprintf("OOB R^2 (PC1)=%.3f\n", oob_r2(rf_pc1, dat$PC1)))
cat(sprintf("OOB R^2 (PC2)=%.3f\n", oob_r2(rf_pc2, dat$PC2)))

make_importance_df <- function(model, axis_name) {
  data.frame(
    var = names(ranger::importance(model)),
    imp = as.numeric(ranger::importance(model)),
    axis = axis_name,
    row.names = NULL
  )
}

imp_all <- rbind(
  make_importance_df(rf_pc1, "PC1"),
  make_importance_df(rf_pc2, "PC2")
)

# Replace negative importance values with zero for visualization
imp_all$imp[imp_all$imp < 0] <- 0

# Fix variable order
imp_all$axis <- factor(imp_all$axis, levels = c("PC1", "PC2"))
imp_all$var <- factor(
  imp_all$var,
  levels = c("TDP", "TDN", "BIX", "HIX", "POC", "DOC", "d13C")
)

# Save importance table
readr::write_csv(
  imp_all,
  paste0("outputs/rf_importance_PC1_PC2_", dataset, ".csv")
)

# ============================================================
# Plot Random Forest Results
# ============================================================
p_rf <- ggplot(imp_all, aes(x = var, y = imp, fill = axis)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_flip() +
  scale_x_discrete(
    labels = c(
      d13C = expression(delta^{13}*C),
      DOC = "DOC",
      POC = "POC",
      TDN = "TDN",
      TDP = "TDP",
      HIX = "HIX",
      BIX = "BIX"
    )
  ) +
  scale_fill_manual(values = c(
    "PC1" = "#E59880",
    "PC2" = "#5A6FAE"
  )) +
  labs(
    x = NULL,
    y = "Permutation importance",
    fill = "Axis"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    aspect.ratio = 1.5,
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.7)
  )

print(p_rf)

ggsave(
    filename = "Fig5A.pdf",
    plot = p_rf,
    device = cairo_pdf,
	width = 110,
    height = 180,
    units = "mm",
    bg = "white"
)
