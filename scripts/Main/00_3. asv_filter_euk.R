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
sample_sheet <- read.csv("data/eukaryote/metadata_euk.csv", row.names = 1)
asv_sheet <- read.csv("data/eukaryote/otu_table.csv", row.names = 1)
tax_sheet <- read.csv("data/eukaryote/taxonomy_table.csv", row.names = 1)

# ============================================================
# Clean taxonomy table
# ============================================================
std_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

missing_cols <- setdiff(std_cols, colnames(tax_sheet))
for (cl in missing_cols) tax_sheet[[cl]] <- NA_character_
tax_sheet <- tax_sheet[, std_cols, drop = FALSE]

strip_prefix <- function(x) {
  x <- as.character(x)
  x <- gsub("^[dkpcofgs]__", "", x, perl = TRUE)
  x <- gsub("^D_[0-9]+__", "", x, perl = TRUE)
  x
}

tax_sheet[] <- lapply(tax_sheet, strip_prefix)

tax_sheet[] <- lapply(tax_sheet, function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- "Unassigned"
  x
})

tax_sheet$Kingdom <- case_when(
  str_detect(str_to_lower(tax_sheet$Kingdom), "eukaryot|eukarya") ~ "Eukaryota",
  str_detect(str_to_lower(tax_sheet$Kingdom), "bacteria") ~ "Bacteria",
  str_detect(str_to_lower(tax_sheet$Kingdom), "archaea") ~ "Archaea",
  str_to_lower(tax_sheet$Kingdom) == "unassigned" ~ "Unassigned",
  TRUE ~ tax_sheet$Kingdom
)

otu_mat <- as.matrix(asv_sheet)

missing_taxa <- setdiff(colnames(otu_mat), rownames(tax_sheet))

if (length(missing_taxa) > 0) {
  tax_sheet[missing_taxa, std_cols] <- "Unassigned"
}

tax_sheet <- tax_sheet[colnames(otu_mat), std_cols, drop = FALSE]

# ============================================================
# Construct phyloseq object
# ============================================================
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# ============================================================
# Remove non-target ASVs
# ============================================================
`%ni%` <- Negate(`%in%`)

drop_animals <- c(
  "Metazoa", "Arthropoda", "Annelida", "Mollusca", "Platyhelminthes",
  "Cnidaria", "Vertebrata", "Rotifera", "Gastrotricha", "Nematozoa",
  "Nematoda", "Tardigrada", "Ctenophora", "Porifera", "Echinodermata",
  "Hemichordata", "Bryozoa", "Brachiopoda", "Onychophora", "Priapulida",
  "Kinorhyncha", "Chaetognatha", "Gnathostomulida", "Loricifera", "Placozoa"
)

drop_plants <- c(
  "Embryophyta", "Streptophyta", "Phragmoplastophyta",
  "Tracheophyta", "Klebsormidiophyceae"
)

drop_fungi <- c(
  "Fungi", "Ascomycota", "Basidiomycota", "Chytridiomycota",
  "Zoopagomycota", "Cryptomycota", "Rozellomycota",
  "Mucoromycota", "Blastocladiomycota", "Protosporangiida",
  "Peronosporomycetes"
)

drop_organelle_class <- c("Chloroplast","Mitochondria")

ps18 <- subset_taxa(
  ps_all,
  Kingdom %in% c("Eukaryota", "Unassigned") &
    (Phylum %ni% c(drop_animals, drop_plants, drop_fungi)) &
    (Class %ni% drop_organelle_class)
)

# ============================================================
# Select sample types
# ============================================================
exp <- subset_samples(ps18, Type %in% c("leaf", "tile_l", "tile_d", "sediment_5", "sediment_20"))

# ============================================================
# Sample quality control
# ============================================================
sc_threshold        <- 0.85
min_reads           <- 1000
min_obs_asv         <- 50

mean_rel_abund_th   <- 0.0005
min_total_count_asv <- 10

stopifnot(inherits(exp, "phyloseq"))

# Calculate sample coverage
otu_raw <- otu_table(exp)
mat_raw <- as(otu_raw, "matrix")
if (!taxa_are_rows(otu_raw)) mat_raw <- t(mat_raw)

reads_raw <- colSums(mat_raw, na.rm = TRUE)
obs_asv   <- colSums(mat_raw > 0, na.rm = TRUE)

sc_per_sample_raw <- sapply(seq_len(ncol(mat_raw)), function(j){
  if (reads_raw[j] <= 0) return(NA_real_)
  as.numeric(iNEXT::DataInfo(mat_raw[, j], datatype = "abundance")$SC)
})
names(sc_per_sample_raw) <- colnames(mat_raw)

# Retain high-quality samples
keep_samples_base <- (!is.na(sc_per_sample_raw)) &
                     (sc_per_sample_raw >= sc_threshold) &
                     (reads_raw >= min_reads) &
                     (obs_asv   >= min_obs_asv)

rescue <- (!is.na(sc_per_sample_raw)) &
          (sc_per_sample_raw >= 0.995) &
          (reads_raw >= ceiling(0.8 * min_reads)) &
          (obs_asv  >= 10)

keep_samples <- keep_samples_base | rescue
keep_samples[is.na(keep_samples)] <- FALSE

kept_ids  <- names(sc_per_sample_raw)[keep_samples %in% TRUE]
exp_kept <- prune_samples(kept_ids, exp)

# Filter low-abundance ASVs
otu2 <- otu_table(exp_kept)
mat2 <- as(otu2, "matrix")
if (!taxa_are_rows(otu2)) mat2 <- t(mat2)

asv_total <- rowSums(mat2, na.rm = TRUE)

colsum2 <- colSums(mat2, na.rm = TRUE); colsum2[colsum2 == 0] <- 1
rel2    <- sweep(mat2, 2, colsum2, "/")
rel2[!is.finite(rel2)] <- 0
mean_rel2_pos <- rowSums(rel2, na.rm = TRUE) / pmax(1, rowSums(mat2 > 0, na.rm = TRUE))

keep_taxa2 <- (asv_total >= min_total_count_asv) &
              (mean_rel2_pos >= mean_rel_abund_th)

ps_final <- prune_taxa(rownames(mat2)[keep_taxa2], exp_kept)

# Generate QC summary
sc_table <- data.frame(
  SampleID = names(sc_per_sample_raw),
  Reads    = as.integer(reads_raw),
  ObsASV   = as.integer(obs_asv),
  SC_raw   = as.numeric(sc_per_sample_raw),
  Keep     = keep_samples
)[order(keep_samples, decreasing = TRUE), ]

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

# ============================================================
# Classify eukaryotic ASVs into algae and protists
# ============================================================
stopifnot(inherits(ps_final, "phyloseq"))

otu_m <- as(otu_table(ps_final), "matrix")
if (!taxa_are_rows(otu_table(ps_final))) otu_m <- t(otu_m)

tx <- as.data.frame(tax_table(ps_final), stringsAsFactors = FALSE)
tx[] <- lapply(tx, function(x){ x <- trimws(as.character(x)); x[x=="" | is.na(x)] <- "Unassigned"; x })

algae_phyla <- c("Diatomea","Bacillariophyta","Chlorophyta","Ochrophyta",
                 "Cryptophyceae","Prymnesiophyceae","Haptophyta","Florideophycidae")
				 
protist_phyla <- c("Ciliophora","Cercozoa","Amoebozoa","Retaria","Centrohelida",
                   "Heterolobosea","Bicosoecida","Euglenozoa","Dinoflagellata",
                   "Apicomplexa")

is_unknown <- function(p){
  p <- tolower(p)
  p %in% tolower(c("Unassigned","Incertae_Sedis","SAR")) |
    grepl("uncultured|unknown|metagenome|environmental|MAST", p, ignore.case=TRUE)
}

tx$AP_Group <- ifelse(tx$Phylum %in% algae_phyla,   "Algae",
                 ifelse(tx$Phylum %in% protist_phyla,"Protists",
                 ifelse(is_unknown(tx$Phylum),       "Unknown", "Other")))

# ============================================================
# Summarize algae/protist read composition
# ============================================================
grp <- rowsum(otu_m, group = tx$AP_Group)

missing_groups <- setdiff(c("Algae", "Protists", "Unknown", "Other"), rownames(grp))
if (length(missing_groups) > 0) {
  grp <- rbind(
    grp,
    matrix(
      0,
      nrow = length(missing_groups),
      ncol = ncol(grp),
      dimnames = list(missing_groups, colnames(grp))
    )
  )
}

grp <- grp[c("Algae", "Protists", "Unknown", "Other"), , drop = FALSE]
reads_tbl <- data.frame(
  SampleID       = colnames(grp),
  Algae_reads    = as.integer(grp["Algae", ]),
  Protists_reads = as.integer(grp["Protists", ]),
  Unknown_reads  = as.integer(grp["Unknown", ]),
  Other_reads    = as.integer(grp["Other", ])
)

reads_tbl$total_reads <- with(
  reads_tbl,
  Algae_reads + Protists_reads + Unknown_reads + Other_reads
)

reads_tbl$prop_Algae <- ifelse(
  reads_tbl$total_reads > 0,
  reads_tbl$Algae_reads / reads_tbl$total_reads,
  0
)

reads_tbl$prop_Protists <- ifelse(
  reads_tbl$total_reads > 0,
  reads_tbl$Protists_reads / reads_tbl$total_reads,
  0
)

reads_tbl$prop_Unknown <- ifelse(
  reads_tbl$total_reads > 0,
  reads_tbl$Unknown_reads / reads_tbl$total_reads,
  0
)

reads_tbl$prop_Other <- 1 - (
  reads_tbl$prop_Algae +
    reads_tbl$prop_Protists +
    reads_tbl$prop_Unknown
)

# ============================================================
# Create algae/protist datasets for downstream analyses
# ============================================================
algae_taxa    <- rownames(tx)[tx$AP_Group == "Algae"]
protists_taxa <- rownames(tx)[tx$AP_Group == "Protists"]
ap_taxa       <- rownames(tx)[tx$AP_Group %in% c("Algae", "Protists")]

ps_algae <- prune_taxa(algae_taxa, ps_final)
ps_algae <- prune_samples(sample_sums(ps_algae) > 0, ps_algae)

ps_protists <- prune_taxa(protists_taxa, ps_final)
ps_protists <- prune_samples(sample_sums(ps_protists) > 0, ps_protists)

ps_ap <- prune_taxa(ap_taxa, ps_final)
ps_ap <- prune_samples(sample_sums(ps_ap) > 0, ps_ap)

saveRDS(exp_kept, "outputs/exp_kept_euk.rds")
saveRDS(ps_final, "outputs/ps_final_euk.rds")
saveRDS(ps_algae, "outputs/ps_final_algae.rds")
saveRDS(ps_protists, "outputs/ps_final_protists.rds")
saveRDS(ps_ap, "outputs/ps_final_ap.rds")