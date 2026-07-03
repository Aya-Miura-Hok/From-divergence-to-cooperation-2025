#!/bin/bash
set -euo pipefail

# ============================================================
# Calculate network topology metrics
# ============================================================
# Required input files in data/network_analysis/outputs/:
#   exported-corr/
#   *_corr_export/
#   all_asv_domain.tsv
#   site*_asv_domain.tsv
#
# Required Python scripts (same directory):
#   net_metrics_from_corr.py
#   net_metrics_domains.py

# ============================================================
# Network metrics (overall)
# ============================================================
cd data/network_analysis

python net_metrics_from_corr.py \
  --corr outputs/exported-corr/pairwise_comparisons.tsv \
  --output outputs/all_metrics_signed.csv \
  --min_r 0.3

# ============================================================
# Network metrics (site-specific)
# ============================================================
for site in site1 site2 site3 site4; do
  python net_metrics_from_corr.py \
    --corr outputs/${site}_corr_export/pairwise_comparisons.tsv \
    --output outputs/${site}_metrics_signed.csv \
    --min_r 0.3
done

# ============================================================
# Network metrics (domain-specific)
# ============================================================
python net_metrics_domains.py \
  --corr outputs/exported-corr/pairwise_comparisons.tsv \
  --asv_domain outputs/all_asv_domain.tsv \
  --output outputs/domain_metrics_all.csv \
  --min_r 0.3

for site in site1 site2 site3 site4; do
  python net_metrics_domains.py \
    --corr outputs/${site}_corr_export/pairwise_comparisons.tsv \
    --asv_domain outputs/${site}_asv_domain.tsv \
    --output outputs/domain_metrics_${site}.csv \
    --min_r 0.3
done