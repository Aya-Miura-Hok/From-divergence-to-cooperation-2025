# ============================================================
# Input
# ============================================================
# This script calculates alpha diversity using domain-specific
# filtered phyloseq objects generated in Script 01_preprocessing.R.
#
# Required input files:
#   outputs/exp_kept_bac.rds
#   outputs/exp_kept_fun.rds
#   outputs/exp_kept_euk.rds

# ============================================================
# Load libraries
# ============================================================
library(tidyverse)
library(phyloseq)
library(iNEXT)
library(rstatix)
library(ggpubr)

# Set default ggplot theme
theme_set(theme_bw())

# ============================================================
# Load phyloseq object
# ============================================================
# Select dataset
dataset <- "bac"
# dataset <- "fun"
# dataset <- "euk"

exp_kept <- readRDS(paste0("outputs/exp_kept_", dataset, ".rds"))

# ============================================================
# Convert each sample to an abundance vector for iNEXT
# ============================================================
otu <- as(otu_table(exp_kept), "matrix"); if(!taxa_are_rows(exp_kept)) otu <- t(otu)
sids <- sample_names(exp_kept)
samp_vecs <- setNames(lapply(sids, function(s){ x <- otu[,s]; x[x>0] }), sids)

meta <- data.frame(sample_data(exp_kept))
meta$SampleID <- rownames(meta)

sc <- iNEXT::DataInfo(samp_vecs, "abundance")$SC
sc_df <- data.frame(SampleID=names(samp_vecs), SC_obs=sc)
meta_sc <- left_join(meta, sc_df, by="SampleID")

sc_star_by_site <- meta_sc %>%
  group_by(Site) %>% summarise(SC_star = pmin(min(SC_obs, na.rm=TRUE), 0.999), .groups="drop")

# ============================================================
# Calculate alpha diversity using site-specific sample coverage
# ============================================================
alpha_siteSC <- map_dfr(unique(meta_sc$Site), function(st){
  lv  <- sc_star_by_site$SC_star[sc_star_by_site$Site==st]
  ids <- meta_sc %>% filter(Site==st, is.finite(SC_obs)) %>% pull(SampleID)
  if(!length(ids)) return(NULL)
  est <- iNEXT::estimateD(samp_vecs[ids], q=c(0,1,2), datatype="abundance",
                          base="coverage", level=lv, nboot=300)
  est$Site <- st; est$SC_star <- lv; est
}) %>%
  rename(SampleID=Assemblage) %>%
  left_join(meta_sc[,c("SampleID","Site","Type")], by=c("SampleID","Site")) %>%
  select(SampleID, Site, Type, Order.q, qD, qD.LCL, qD.UCL, SC_star)

# ============================================================
# Save alpha diversity results
# ============================================================
outdir <- "results_inext"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

write.csv(
  alpha_siteSC,
  file.path(outdir, paste0("alpha_siteSC_", dataset, "_long.csv")),
  row.names = FALSE
)

csv_file <- file.path(outdir, paste0("alpha_siteSC_", dataset, "_long.csv"))

# ============================================================
# Prepare data for statistical analysis and plotting
# ============================================================
df <- read_csv(csv_file, show_col_types = FALSE)

# Plot Shannon diversity (q = 1)
q_show <- 1
dfq <- df %>% filter(Order.q == q_show)

# Keep only sites containing at least two substrate types
valid_sites <- dfq %>% distinct(Site, Type) %>%
  count(Site, name = "n_types") %>%
  filter(n_types >= 2) %>% pull(Site)

dfq <- dfq %>% filter(Site %in% valid_sites)

# ============================================================
# Statistical analysis
# ============================================================
dunn_tbl <- dfq %>%
  dplyr::group_by(Site) %>%
  rstatix::dunn_test(qD ~ Type, p.adjust.method = "BH") %>%
  rstatix::add_xy_position(x = "Type") %>%
  dplyr::ungroup()

dunn_sig <- dunn_tbl %>%
  dplyr::filter(p.adj < 0.05)

# ============================================================
# Plot alpha diversity
# ============================================================
dunn_sig_plot <- dunn_sig %>%
  ungroup()

p <- ggplot(dfq, aes(x = Type, y = qD, fill = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  facet_wrap(
    ~ Site,
    scales = "free_y",
    labeller = labeller(
      Site = c(
        "1" = "Site 1",
        "2" = "Site 2",
        "3" = "Site 3",
        "4" = "Site 4"
      )
    )
  ) +
  { if (nrow(dunn_sig_plot) > 0)
      ggpubr::stat_pvalue_manual(
        dunn_sig_plot,
        label       = "p.adj.signif",
        y.position  = "y.position",
        xmin        = "xmin",
        xmax        = "xmax",
        size        = 6,
        inherit.aes = FALSE
      )
    else NULL } +
  labs(
    x = "Substrate type",
    y = "Shannon diversity"
  ) +
  theme_minimal(base_size = 15, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      size = 16,
      face = "bold"
    )
  )

print(p)