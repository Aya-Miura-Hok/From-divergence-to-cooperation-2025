# ============================================================
# Input
# ============================================================
# This script uses exp_kept_* objects 
# generated in Script 00_1. asv_filter_bac.R, 00_2. asv_filter_fun.R, 00_3. asv_filter_euk.R for gamma diversity analyses.
# These objects contain samples that passed sample-level QC before ASV-level filtering.
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

theme_set(theme_bw())

# ============================================================
# Calculate γ-diversity by site
# ============================================================
pool_by_site <- function(ps) {
  sd <- data.frame(sample_data(ps))
  sd$SampleID <- rownames(sd)

  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  sp <- split(sd$SampleID, sd$Site)

  setNames(
    lapply(sp, function(v) {
      rs <- rowSums(otu[, v, drop = FALSE])
      rs[rs > 0]
    }),
    names(sp)
  )
}

calc_gamma_site <- function(ps) {
  glist_site <- pool_by_site(ps)

  sc_site <- iNEXT::DataInfo(glist_site, datatype = "abundance")$SC
  sc_star <- pmin(min(sc_site, na.rm = TRUE), 0.999)

  iNEXT::estimateD(
    glist_site,
    q = c(0, 1, 2),
    datatype = "abundance",
    base = "coverage",
    level = sc_star,
    nboot = 300,
    conf = 0.95
  ) %>%
    rename(Site = Assemblage)
}

# Load domain-specific filtered phyloseq objects
exp_kept_bac <- readRDS("outputs/exp_kept_bac.rds")
exp_kept_fun <- readRDS("outputs/exp_kept_fun.rds")
exp_kept_euk <- readRDS("outputs/exp_kept_euk.rds")

# Calculate gamma diversity
bac_gamma_site <- calc_gamma_site(exp_kept_bac)
fun_gamma_site <- calc_gamma_site(exp_kept_fun)
euk_gamma_site <- calc_gamma_site(exp_kept_euk)

# Save results
saveRDS(bac_gamma_site, file = "outputs/bac_gamma_site.rds")
saveRDS(fun_gamma_site, file = "outputs/fun_gamma_site.rds")
saveRDS(euk_gamma_site, file = "outputs/euk_gamma_site.rds")

# ============================================================
# plot results
# ============================================================
gamma_all <- bind_rows(
  bac_gamma_site %>% mutate(Domain = "Bacteria"),
  fun_gamma_site %>% mutate(Domain = "Fungi"),
  euk_gamma_site %>% mutate(Domain = "Eukaryote")
)

gamma_all$Domain <- factor(
  gamma_all$Domain,
  levels = c("Bacteria", "Fungi", "Eukaryote")
)

p_gamma <- ggplot(
  gamma_all,
  aes(x = Site, y = qD, color = Domain, group = Domain)
) +
  geom_pointrange(
    aes(ymin = qD.LCL, ymax = qD.UCL),
    position = position_dodge(width = 0.4),
    size = 0.7
  ) +
  facet_wrap(
    ~Order.q,
    scales = "free_y",
    labeller = labeller(
      Order.q = c("0" = "q = 0", "1" = "q = 1", "2" = "q = 2")
    )
  ) +
  scale_color_manual(
    name = "Microbial group",
    values = c(
      "Bacteria" = "#E64B35FF",
      "Fungi" = "#4DBBD5FF",
      "Eukaryote" = "#00A087FF"
    )
  ) +
  labs(
    x = "Site",
    y = expression(q[D]~"(estimate ± 95% CI)")
  ) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

print(p_gamma)