#!/usr/bin/env python
import argparse, pandas as pd, numpy as np, networkx as nx

def normalize_domain(x):
    if pd.isna(x): return "Unassigned"
    x = str(x)
    x = x.replace("d__","").replace("k__","")
    x = x.strip()
    # 好みでまとめる
    if x.lower().startswith("bact"): return "Bacteria"
    if x.lower().startswith("fung"): return "Fungi"
    if x.lower().startswith("euk"):  return "Eukaryota"
    if x.lower().startswith("prot"): return "Protists"
    if x.lower().startswith("alga"): return "Algae"
    return x

def metrics_from_graph(G):
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    density = nx.density(G) if n_nodes>1 else 0.0
    trans   = nx.transitivity(G) if n_edges>0 else np.nan
    mean_deg= np.mean([d for _,d in G.degree()]) if n_nodes>0 else 0.0
    # 正負比率（全エッジ対象）
    ws = [d.get("weight",0.0) for _,_,d in G.edges(data=True)]
    if len(ws)>0:
        pos_ratio = np.mean([w>0 for w in ws])
        neg_ratio = 1.0 - pos_ratio
    else:
        pos_ratio = neg_ratio = np.nan
    # モジュラリティは正エッジのみで（community-louvainが負に非対応のため）
    try:
        import community as community_louvain
        Gpos = nx.Graph([(u,v,d) for u,v,d in G.edges(data=True) if d.get("weight",0)>0])
        if Gpos.number_of_edges()>0:
            part = community_louvain.best_partition(Gpos, weight="weight")
            modularity = community_louvain.modularity(part, Gpos, weight="weight")
        else:
            modularity = np.nan
    except Exception:
        modularity = np.nan
    # ハブ割合（上位5%）
    degs = np.array([d for _,d in G.degree()])
    if len(degs)>0:
        cutoff = np.percentile(degs, 95)
        hub_rel5p = np.mean(degs >= cutoff)
    else:
        hub_rel5p = np.nan
    return dict(n_nodes=n_nodes,n_edges=n_edges,density=density,
                transitivity=trans,mean_deg=mean_deg,
                pos_ratio=pos_ratio,neg_ratio=neg_ratio,
                modularity=modularity,hub_rel5p=hub_rel5p)

def build_subset_graph(pairs, domap, left_set, right_set, min_r):
    # 端点ドメイン条件を満たすエッジのみ採用
    rows=[]
    for _,r in pairs.iterrows():
        s,t,rv = r["Source"], r["Target"], r["R"]
        if abs(rv) < min_r: continue
        ds = domap.get(s, "Unassigned")
        dt = domap.get(t, "Unassigned")
        ok = ((ds in left_set and dt in right_set) or
              (ds in right_set and dt in left_set))
        if ok:
            rows.append((s,t,rv))
    G = nx.Graph()
    for s,t,rv in rows:
        G.add_edge(s,t,weight=float(rv))
    return G

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--corr", required=True, help="pairwise_comparisons.tsv")
    ap.add_argument("--asv_domain", required=True, help="ASV→Domain tsv (ASV,Kingdom)")
    ap.add_argument("--output", required=True, help="output csv")
    ap.add_argument("--min_r", type=float, default=0.3, help="abs(r) cutoff")
    # ドメイン名のグルーピング（必要に応じて調整）
    args = ap.parse_args()

    # 相関ペア
    df = pd.read_csv(args.corr, sep="\t")
    df = df.rename(columns={df.columns[0]:"Source", df.columns[1]:"Target", df.columns[2]:"R"})

    # ASV→ドメイン
    dm = pd.read_csv(args.asv_domain, sep="\t")
    # 列名の揺れ吸収
    if "ASV" not in dm.columns:
        dm.rename(columns={dm.columns[0]:"ASV"}, inplace=True)
    if "Kingdom" not in dm.columns:
        dm.rename(columns={dm.columns[1]:"Kingdom"}, inplace=True)
    dm["DomainN"] = dm["Kingdom"].map(normalize_domain)
    domap = dict(zip(dm["ASV"], dm["DomainN"]))

    # カテゴリ定義（調整可）
    domCats = {
        "Bac–Bac": ({"Bacteria"}, {"Bacteria"}),
        "Fun–Fun": ({"Fungi"}, {"Fungi"}),
        "Bac–Fun": ({"Bacteria"}, {"Fungi"}),
        "Bac–Euk": ({"Bacteria"}, {"Eukaryota","Protists","Algae"}),
        "Fun–Euk": ({"Fungi"}, {"Eukaryota","Protists","Algae"}),
        "Euk–Euk": ({"Eukaryota","Protists","Algae"}, {"Eukaryota","Protists","Algae"}),
    }

    out_rows=[]
    for label, (L,R) in domCats.items():
        G = build_subset_graph(df, domap, L, R, args.min_r)
        m = metrics_from_graph(G)
        m.update({"category": label})
        out_rows.append(m)

    pd.DataFrame(out_rows, columns=["category","n_nodes","n_edges","density",
                                    "transitivity","mean_deg","pos_ratio","neg_ratio",
                                    "modularity","hub_rel5p"]).to_csv(args.output, index=False)
    print("✅ saved:", args.output)