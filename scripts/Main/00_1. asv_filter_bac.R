# ============================================================
# Load libraries
# ============================================================
library(tidyverse)
library(phyloseq)
library(iNEXT)

# Set default ggplot theme
theme_set(theme_bw())

# ============================================================
# Load input data (metadata, OTU table, taxonomy)
# ============================================================
sample_sheet <- read.csv("data/bacteria/metadata_bac.csv", row.names = 1)
asv_sheet <- read.csv("data/bacteria/otu_table.csv", row.names = 1)
asv_sheet <- t(asv_sheet)
tax_sheet <- read.csv("data/bacteria/taxonomy_table.csv", row.names = 1)

# ============================================================
# Construct phyloseq object
# ============================================================
ps_all <- phyloseq(otu_table(asv_sheet, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))

ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")

# ============================================================
# Remove non-target ASVs
# ============================================================
tax <- tax_table(ps)
tax_mat <- as.matrix(tax)

# Remove chloroplast, mitochondrial, and non-target sequences
is_organelle <- apply(tax_mat, 1, function(x) {
  any(grepl("chloroplast|mitochondria", x, ignore.case = TRUE), na.rm = TRUE)
})
is_euk <- FALSE
if ("Kingdom" %in% colnames(tax_mat)) {
  is_euk <- tax_mat[, "Kingdom"] %in% c("Eukaryota","d__Eukaryota")
}

drop_taxa <- taxa_names(ps)[is_organelle | is_euk]
ps_no_plastid <- prune_taxa(setdiff(taxa_names(ps), drop_taxa), ps)
exp <- subset_samples(ps_no_plastid, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))

# ============================================================
# Sample quality control
# ============================================================
# Filtering parameters
sc_threshold       <- 0.90
min_reads          <- 1000
min_obs_asv        <- 50
prev_prop          <- 0.10
mean_rel_abund_th  <- 0.0005
min_total_count_asv <- 10

stopifnot(inherits(exp, "phyloseq"))

# Calculate sample coverage
otu_raw <- otu_table(exp)
mat_raw <- as(otu_raw, "matrix")
if (!taxa_are_rows(otu_raw)) mat_raw <- t(mat_raw) # rows = ASVs, columns = samples

reads_raw <- colSums(mat_raw, na.rm = TRUE)
obs_asv   <- colSums(mat_raw > 0, na.rm = TRUE)

sc_per_sample_raw <- sapply(1:ncol(mat_raw), function(j){
  if (reads_raw[j] <= 0) return(NA_real_)
  as.numeric(DataInfo(mat_raw[, j], datatype = "abundance")$SC)
})
names(sc_per_sample_raw) <- colnames(mat_raw)

# Retain high-quality samples
keep_samples <- (!is.na(sc_per_sample_raw)) &
                (sc_per_sample_raw >= sc_threshold) &
                (reads_raw >= min_reads) &
                (obs_asv   >= min_obs_asv)

exp_kept <- prune_samples(names(keep_samples)[keep_samples], exp)

# Filter low-abundance ASVs
otu2 <- otu_table(exp_kept)
mat2 <- as(otu2, "matrix")
if (!taxa_are_rows(otu2)) mat2 <- t(mat2)

nsamp2 <- ncol(mat2)

asv_total <- rowSums(mat2, na.rm = TRUE)

prev2 <- rowSums(mat2 > 0, na.rm = TRUE) / nsamp2
rel2        <- sweep(mat2, 2, colSums(mat2, na.rm = TRUE), "/")
rel2[!is.finite(rel2)] <- 0
mean_rel2   <- rowMeans(rel2, na.rm = TRUE)

keep_taxa2 <- (asv_total >= min_total_count_asv) &
              (prev2 >= prev_prop) &
              (mean_rel2 >= mean_rel_abund_th)

ps_final <- prune_taxa(rownames(mat2)[keep_taxa2], exp_kept)

# Generate QC summary
sc_table <- data.frame(
  SampleID = names(sc_per_sample_raw),
  Reads    = as.integer(reads_raw),
  ObsASV   = as.integer(obs_asv),
  SC_raw   = as.numeric(sc_per_sample_raw),
  Keep     = keep_samples
)
sc_table <- sc_table[order(sc_table$Keep, decreasing = TRUE), ]

# Check replicate numbers
meta <- data.frame(sample_data(ps_final))
group_cols <- c("Site", "Type")

rep_count <- meta %>%
  group_by(across(all_of(group_cols))) %>%
  summarise(n_rep = n(), .groups = "drop")

print(rep_count)

if(any(rep_count$n_rep < 2)){
  message("Some groups have fewer than two replicates.")
  print(filter(rep_count, n_rep < 2))
} else {
  message("All groups have at least two replicates.")
}

# Create datasets for downstream analyses
saveRDS(exp_kept, "outputs/exp_kept_bac.rds")
saveRDS(ps_final, "outputs/ps_final_bac.rds")
ps_alpha <- exp_kept
ps_beta  <- ps_final
