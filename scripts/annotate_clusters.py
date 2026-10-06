import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mine = pd.read_csv("data/obs_neurons.csv", index_col=0)
rds = pd.read_csv("data/neurons_clustered_metadata.csv")
rds_genes = set(open("data/rds_genes.txt").read().split())

rows = []
for sp_, tag in [("P.man", "Pman"), ("P.pol", "Ppol")]:
    a = sc.read_h5ad(f"data/{tag}_raw.h5ad")
    a = a[a.obs_names.isin(mine.index)].copy()
    nm = a.var["gene_name"].astype(str)
    keep = nm.isin(rds_genes).values & ~nm.duplicated(keep="first").values
    X = a.X[:, np.where(keep)[0]]
    rows.append(pd.DataFrame({"nCount": np.asarray(X.sum(1)).ravel(),
                              "nFeature": np.asarray((X > 0).sum(1)).ravel()}, index=a.obs_names))
mine = mine.join(pd.concat(rows))
mine["sp"] = mine["species"].map({"P.man": "BW", "P.pol": "PO"})
mine["key"] = mine["sp"] + "|" + mine["sex"] + "|" + mine["nCount"].astype(int).astype(str) + "|" + mine["nFeature"].astype(int).astype(str)
rds["key"] = rds["species"] + "|" + rds["sex"] + "|" + rds["nCount_RNA"].astype(int).astype(str) + "|" + rds["nFeature_RNA"].astype(int).astype(str)
u_m = mine["key"].map(mine["key"].value_counts()) == 1
u_r = rds["key"].map(rds["key"].value_counts()) == 1
m = mine[u_m].reset_index().merge(rds[u_r][["key", "rep", "cluster", "UMAP_1", "UMAP_2"]],
                                  on="key", suffixes=("", "_rds"))
print(f"mine {len(mine)}, authors {len(rds)}, matched one-to-one {len(m)} ({len(m)/len(rds):.1%} of authors' neurons)")
print(pd.crosstab(m["rep"], m["rep_rds"]))
m.to_csv("data/matched_nuclei.csv", index=False)

a_, b_ = m["cluster"].astype(str), m["cluster_rds"].astype(str)
print("ARI:", round(ARI(a_, b_), 3), " NMI:", round(NMI(a_, b_), 3))
ct = pd.crosstab(a_, b_)
ct.to_csv("results/contingency_vs_authors.csv")
purity = ct.max(axis=1) / ct.sum(axis=1)
print("my clusters with >=50% of nuclei in one authors' cluster:", int((purity >= .5).sum()), "of", len(purity))
print("median purity:", round(purity.median(), 2))

rn = ct.div(ct.sum(1), axis=0)
order_r = rn.idxmax(axis=1).astype(int).sort_values().index
rn = rn.loc[order_r, sorted(rn.columns, key=int)]
plt.figure(figsize=(12, 10))
sns.heatmap(rn, cmap="viridis", xticklabels=True, yticklabels=True, cbar_kws={"label": "fraction of my cluster"})
plt.xlabel("Authors' cluster (RDS)"); plt.ylabel("My cluster")
plt.tight_layout(); plt.savefig("figures/supplemental/contingency_vs_authors.png", dpi=130)

np.random.seed(0)
MAX_PER_CLUSTER = 3000
FRACTION_RULE = 0.15

ref_all = pd.read_csv("data/moffitt_merfish.csv", dtype={"Neuron_cluster_ID": "string"})
print("MERFISH cells:", len(ref_all), "| behaviours:", ref_all["Behavior"].value_counts().to_dict())
ref_all = ref_all[ref_all["Neuron_cluster_ID"].notna() & (ref_all["Cell_class"] != "Ambiguous")]
ref_all = ref_all[ref_all["Behavior"] == "Naive"]
meta_cols = ["Cell_ID", "Animal_ID", "Animal_sex", "Behavior", "Bregma", "Centroid_X", "Centroid_Y",
             "Cell_class", "Neuron_cluster_ID"]
gene_cols = [c for c in ref_all.columns if c not in meta_cols and not c.startswith("Blank")]
gene_cols = [c for c in gene_cols if not ref_all[c].isna().any()]
print("naive neurons:", len(ref_all), "| reference genes:", len(gene_cols))
ref_all["type"] = np.where(ref_all["Neuron_cluster_ID"].str.startswith("I"), "Inhibitory", "Excitatory")

neu = sc.read_h5ad("data/neurons_integrated.h5ad")
q_all = neu.raw.to_adata()
q_all.obs = neu.obs.copy()
genes = [g for g in gene_cols if g in q_all.var_names]
print("genes shared with my data:", len(genes), "of", len(gene_cols))
missing = sorted(set(gene_cols) - set(genes)); print("not found in my data:", missing)


def scaled(X):
    X = np.asarray(X, dtype=np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-6)
    return np.clip(X, -10, 10)


obs_out, summ = [], []
for ntype in ["Inhibitory", "Excitatory"]:
    r = ref_all[ref_all["type"] == ntype]
    r = r.groupby("Neuron_cluster_ID", group_keys=False).apply(
        lambda x: x.sample(min(len(x), MAX_PER_CLUSTER), random_state=0))
    ref = ad.AnnData(scaled(np.log1p(r[genes].values)), obs=pd.DataFrame(
        {"label": r["Neuron_cluster_ID"].astype(str).values}, index=r["Cell_ID"].astype(str).values))
    ref.var_names = genes
    sc.pp.pca(ref, n_comps=30, random_state=0)
    sc.pp.neighbors(ref, n_neighbors=15, random_state=0)
    sc.tl.umap(ref, random_state=0)

    q = q_all[q_all.obs["neuron_type"] == ntype]
    qa = ad.AnnData(scaled(q[:, genes].X.toarray()), obs=q.obs[["cluster", "species", "rep", "sex"]].copy())
    qa.var_names = genes
    sc.tl.ingest(qa, ref, obs="label", embedding_method="pca", labeling_method="knn")
    qa.obs["moffitt"] = qa.obs["label"].astype(str)

    half = np.random.RandomState(0).rand(ref.n_obs) < 0.5
    a, b = ref[half].copy(), ref[~half].copy()
    sc.pp.pca(a, n_comps=30, random_state=0); sc.pp.neighbors(a, n_neighbors=15, random_state=0)
    sc.tl.ingest(b, a, obs="label", embedding_method="pca", labeling_method="knn")
    b_true = ref[~half].obs["label"].astype(str).values
    acc = (b.obs["label"].astype(str).values == b_true).mean()
    print(f"{ntype}: reference self-transfer accuracy (held-out half): {acc:.2f}")
    summ.append({"type": ntype, "ref_cells": ref.n_obs, "query_nuclei": qa.n_obs,
                 "ref_clusters": ref.obs['label'].nunique(), "ref_self_accuracy": round(acc, 3)})
    obs_out.append(qa.obs[["cluster", "moffitt", "species"]].assign(neuron_type=ntype))

obs = pd.concat(obs_out).loc[neu.obs_names]
obs.to_csv("data/moffitt_labels_per_nucleus.csv")
pd.DataFrame(summ).to_csv("results/moffitt_transfer_summary.csv", index=False)

frac = pd.crosstab(obs["cluster"], obs["moffitt"], normalize="index")
frac.to_csv("results/moffitt_cluster_fractions.csv")
rows = []
for cl, f in frac.iterrows():
    keep = f[f >= FRACTION_RULE].sort_values(ascending=False)
    rows.append({"cluster": cl, "neuron_type": obs.loc[obs.cluster == cl, "neuron_type"].iloc[0],
                 "n_nuclei": int((obs.cluster == cl).sum()),
                 "moffitt_labels": ";".join(keep.index), "fractions": ";".join(f"{v:.2f}" for v in keep.values),
                 "top_label": f.idxmax(), "top_fraction": round(f.max(), 2)})
tab = pd.DataFrame(rows).sort_values("cluster", key=lambda s: s.astype(int))
tab.to_csv("results/moffitt_cluster_labels_15pct.csv", index=False)
print(tab.to_string(index=False))
print("clusters with >=1 label at >=15%:", int((tab.moffitt_labels != "").sum()), "of", len(tab))
print("clusters with exactly 1 label:", int((tab.moffitt_labels.str.count(";") == 0)[tab.moffitt_labels != ""].sum()))
