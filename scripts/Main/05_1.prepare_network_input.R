# ============================================================
# Input
# ============================================================
# Required input files in data/network_analysis/data/:
# bacteria_otu_table.csv, bacteria_taxonomy_table.csv, bacteria_sample_data.csv
# fungi_otu_table.csv, fungi_taxonomy_table.csv, fungi_sample_data.csv
# eukaryote_otu_table.csv, eukaryote_taxonomy_table.csv, eukaryote_sample_data.csv

# ============================================================
# Load libraries
# ============================================================
library(phyloseq)
library(tidyverse)
library(biomformat)

# ============================================================
# Prepare Multi-domain Phyloseq Object for Network Analysis
# ============================================================
# Define Helper Function for Sample Name Cleaning
clean_name <- function(x) {
  x %>%
    str_replace_all("\u3000", " ") %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\.+", "_") %>%
    str_replace_all("-", "_")
}

# Define Helper Function to Reconstruct Domain-specific Phyloseq Objects
read_domain_ps <- function(prefix, domain, dir = "data/network_analysis/data") {
  tax <- read_csv(file.path(dir, paste0(prefix, "_taxonomy_table.csv")), show_col_types = FALSE)
  otu <- read_csv(file.path(dir, paste0(prefix, "_otu_table.csv")),      show_col_types = FALSE)
  sam <- read_csv(file.path(dir, paste0(prefix, "_sample_data.csv")),    show_col_types = FALSE)

  if(!("SampleID" %in% names(sam)) || !("SampleName" %in% names(sam))) {
    stop(sprintf("[%s] sample_data must contain SampleID and SampleName columns.", prefix))
  }
  sam <- sam %>%
    mutate(SampleID = as.character(SampleID),
           SampleName_clean = clean_name(SampleName))
  id2name <- setNames(sam$SampleName_clean, sam$SampleID)

  if(!("ASV" %in% names(otu))) stop(sprintf("[%s] otu_table must contain an ASV column.", prefix))
  otu_df <- column_to_rownames(otu, "ASV")
  colnames(otu_df) <- id2name[colnames(otu_df)]
  keep_cols <- !is.na(colnames(otu_df))
  if(!all(keep_cols)){
    message(sprintf("[%s] Removed unmapped sample columns: %d", prefix, sum(!keep_cols)))
    otu_df <- otu_df[, keep_cols, drop = FALSE]
  }
  otu_mat <- as.matrix(otu_df)
  storage.mode(otu_mat) <- "numeric"

  if(!("ASV" %in% names(tax))) stop(sprintf("[%s] taxonomy_table must contain an ASV column.", prefix))
  tax_df <- column_to_rownames(tax, "ASV")
  tax_df$Domain <- domain
  tax_mat <- as.matrix(tax_df)

  asv_prefix <- paste0(tolower(domain), "_")
  rownames(otu_mat) <- paste0(asv_prefix, rownames(otu_mat))
  rownames(tax_mat) <- paste0(asv_prefix, rownames(tax_mat))

  otu_only <- setdiff(rownames(otu_mat), rownames(tax_mat))
  if(length(otu_only)){
    add_mat <- matrix(NA, nrow=length(otu_only), ncol=ncol(tax_mat),
                      dimnames = list(otu_only, colnames(tax_mat)))
    tax_mat <- rbind(tax_mat, add_mat)
  }
  tax_only <- setdiff(rownames(tax_mat), rownames(otu_mat))
  if(length(tax_only)){
    tax_mat <- tax_mat[setdiff(rownames(tax_mat), tax_only), , drop = FALSE]
  }

  sam_df <- sam %>%
    select(SampleName_clean, everything()) %>%
    column_to_rownames("SampleName_clean")

  common_samples <- intersect(colnames(otu_mat), rownames(sam_df))
  if(length(common_samples) == 0) {
    stop(sprintf("[%s] No shared SampleName was found. Please check sample name formatting.", prefix))
  }
  otu_mat <- otu_mat[, common_samples, drop = FALSE]
  sam_df  <- sam_df[common_samples, , drop = FALSE]

  phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat[rownames(otu_mat), , drop = FALSE]),
    sample_data(sam_df)
  )
}

# Load Domain-specific Datasets
ps_bac <- read_domain_ps(prefix = "bacteria",  domain = "bac")
ps_fun <- read_domain_ps(prefix = "fungi",     domain = "fun")
ps_euk <- read_domain_ps(prefix = "eukaryote", domain = "euk")

# Retain Shared Samples Across Domains
common_samps <- Reduce(intersect, list(sample_names(ps_bac), sample_names(ps_fun), sample_names(ps_euk)))
if(length(common_samps) == 0){
  stop("No shared SampleName was found across the three domains. Please check clean_name().")
}

# Merge Domain-specific Phyloseq Objects
ps_bac_c <- prune_samples(common_samps, ps_bac)
ps_fun_c <- prune_samples(common_samps, ps_fun)
ps_euk_c <- prune_samples(common_samps, ps_euk)

ps_all <- merge_phyloseq(ps_bac_c, ps_fun_c, ps_euk_c)

if (!dir.exists("outputs")) dir.create("outputs")
saveRDS(ps_all, "outputs/ps_all.rds")

# Check Sample Numbers
table(sample_data(ps_all)$Site)

# ============================================================
# Filter ASVs
# ============================================================
otu_mat <- as(otu_table(ps_all), "matrix")
tot     <- rowSums(otu_mat)
prev    <- rowSums(otu_mat > 0) / ncol(otu_mat)

keep    <- (tot >= 10) & (prev >= 0.05)
ps_asv  <- prune_taxa(rownames(otu_mat)[keep], ps_all)

cat(sprintf("ASV kept: %d / %d (%.1f%%)\n", ntaxa(ps_asv), ntaxa(ps_all),
            100*ntaxa(ps_asv)/ntaxa(ps_all)))

# ============================================================
# Export SCNIC Input Files (overall)
# ============================================================
# (A) OTU table → BIOM
otu_mat_out <- as(otu_table(ps_asv), "matrix")
biom_obj    <- make_biom(otu_mat_out)
write_biom(biom_obj, "outputs/otu_table_ASV.biom")

# (B) taxonomy.tsv
tx <- as.data.frame(as(tax_table(ps_asv), "matrix"), stringsAsFactors = FALSE)
tx <- rownames_to_column(tx, var = "Feature ID")

safe <- function(x) ifelse(is.na(x) | x=="", "NA", x)
pref <- function(x, p) paste0(p, safe(x))

taxcols <- colnames(tx)
tx$Taxon <- paste(
  safe(tx[[taxcols[2]]]) |> pref("p__"),
  safe(tx[[taxcols[3]]]) |> pref("c__"),
  safe(tx[[taxcols[4]]]) |> pref("o__"),
  safe(tx[[taxcols[5]]]) |> pref("f__"),
  safe(tx[[taxcols[6]]]) |> pref("g__"),
  safe(tx[[taxcols[7]]]) |> pref("s__"),
  sep="; "
)

write.table(tx[, c("Feature ID","Taxon")],
            file = "outputs/taxonomy_ASV.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# (C) metadata.tsv
md <- as(sample_data(ps_asv), "data.frame")
md <- rownames_to_column(md, var = "RowID")

write.table(md, file = "outputs/metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Exported: otu_table_ASV.biom, taxonomy_ASV.tsv, metadata.tsv\n")

# ============================================================
# Export SCNIC Input Files (site-specific)
# ============================================================
export_ps <- function(ps, prefix){
  stopifnot(inherits(ps, "phyloseq"))
  # OTU (ASV rows x Sample cols)
  otu <- as(otu_table(ps), "matrix"); if(!taxa_are_rows(ps)) otu <- t(otu)
  write_biom(make_biom(otu), paste0(prefix, "_otu.biom"))

  # taxonomy（Feature ID, Taxon）
  tax <- as(tax_table(ps), "matrix") |> as.data.frame()
  tax <- rownames_to_column(tax, "Feature ID")
  safe <- \(x) ifelse(is.na(x)|x=="","NA",x); pref <- \(x,p) paste0(p, safe(x))
  rn <- colnames(tax) # Kingdom Phylum Class Order Family Genus Species
  tax$Taxon <- paste(pref(tax[[rn[2]]],"p__"),pref(tax[[rn[3]]],"c__"),
                     pref(tax[[rn[4]]],"o__"),pref(tax[[rn[5]]],"f__"),
                     pref(tax[[rn[6]]],"g__"),pref(tax[[rn[7]]],"s__"), sep="; ")
  write.table(tax[,c("Feature ID","Taxon")], paste0(prefix, "_taxonomy.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  # metadata
  md <- as(sample_data(ps), "data.frame")
  if(!"SampleID" %in% colnames(md)) md <- rownames_to_column(md, "SampleID")
  write.table(md, paste0(prefix, "_metadata.tsv"), sep="\t",
              row.names=FALSE, quote=FALSE)

  message("exported: ", prefix)
}

ps_site1 <- subset_samples(ps_asv, Site=="1")
ps_site2 <- subset_samples(ps_asv, Site=="2")
ps_site3 <- subset_samples(ps_asv, Site=="3")
ps_site4 <- subset_samples(ps_asv, Site=="4")

export_ps(ps_site1, "outputs/site1")
export_ps(ps_site2, "outputs/site2")
export_ps(ps_site3, "outputs/site3")
export_ps(ps_site4, "outputs/site4")

# ============================================================
# Export ASV-Domain Mapping Tables (overall and site-specific)
# ============================================================
tax <- as.data.frame(as(tax_table(ps_asv), "matrix"))
asv_domain <- tibble::rownames_to_column(tax[, "Kingdom", drop=FALSE], "ASV")
write.table(asv_domain, "outputs/all_asv_domain.tsv", sep="\t", row.names=FALSE, quote=FALSE)

for(s in c("1","2","3","4")){
  ps_s <- get(paste0("ps_site", s))
  tax_s <- as.data.frame(as(tax_table(ps_s), "matrix"))
  asv_domain_s <- tibble::rownames_to_column(tax_s[, "Kingdom", drop=FALSE], "ASV")
  write.table(asv_domain_s, paste0("outputs/site", s, "_asv_domain.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}