#!/usr/bin/env python

import argparse, pandas as pd, numpy as np, networkx as nx
import community as community_louvain

parser = argparse.ArgumentParser()
parser.add_argument("--corr", required=True, help="pairwise_comparisons.tsv")
parser.add_argument("--output", required=True, help="output csv file")
parser.add_argument("--min_r", type=float, default=0.3, help="abs(correlation) cutoff")
args = parser.parse_args()

# --- データ読み込み ---
df = pd.read_csv(args.corr, sep="\t")

# SparCC 出力には 'OTU1','OTU2','r' 列がある前提
df = df.rename(columns={df.columns[0]:"OTU1", df.columns[1]:"OTU2", df.columns[2]:"r"})

# 閾値でフィルタ
df = df[np.abs(df["r"]) >= args.min_r]

# グラフ構築
G = nx.Graph()
for _, row in df.iterrows():
    G.add_edge(row["OTU1"], row["OTU2"], weight=row["r"])

n_nodes = G.number_of_nodes()
n_edges = G.number_of_edges()
density = nx.density(G)
transitivity = nx.transitivity(G)
mean_deg = np.mean([d for _, d in G.degree()]) if n_nodes > 0 else 0

# --- 正負エッジの割合 ---
pos_edges = sum(1 for _,_,d in G.edges(data=True) if d["weight"] > 0)
neg_edges = sum(1 for _,_,d in G.edges(data=True) if d["weight"] < 0)
total_edges = pos_edges + neg_edges if (pos_edges+neg_edges) > 0 else 1

pos_ratio = pos_edges / total_edges
neg_ratio = neg_edges / total_edges

# --- モジュラリティ（正のエッジのみで計算）---
G_pos = nx.Graph([(u,v,d) for u,v,d in G.edges(data=True) if d["weight"] > 0])
if G_pos.number_of_edges() > 0:
    part = community_louvain.best_partition(G_pos, weight="weight")
    modularity = community_louvain.modularity(part, G_pos, weight="weight")
else:
    modularity = np.nan

# --- ハブ相対度数 (上位5%) ---
deg_dict = dict(G.degree())
if len(deg_dict) > 0:
    deg_values = np.array(list(deg_dict.values()))
    cutoff = np.percentile(deg_values, 95)
    hub_rel5p = sum(deg_values >= cutoff) / n_nodes
else:
    hub_rel5p = np.nan

out = pd.DataFrame([{
    "n_nodes": n_nodes,
    "n_edges": n_edges,
    "density": density,
    "transitivity": transitivity,
    "mean_deg": mean_deg,
    "pos_ratio": pos_ratio,
    "neg_ratio": neg_ratio,
    "modularity": modularity,
    "hub_rel5p": hub_rel5p
}])

out.to_csv(args.output, index=False)
print("✅ Network metrics saved to", args.output)