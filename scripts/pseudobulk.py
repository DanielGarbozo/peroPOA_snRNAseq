import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import scanpy as sc

obs = pd.read_csv("data/obs_neurons.csv", index_col=0)
obs["animal"] = obs["species"].str.replace(".", "", regex=False) + "_" + obs["sex"] + "_" + obs["rep"]

ct = pd.crosstab(obs["cluster"], obs["animal"])
ct.to_csv("results/cluster_by_animal_counts.csv")
meta = (obs.drop_duplicates("animal")[["animal", "species", "sex", "rep"]]
        .set_index("animal").loc[ct.columns])
meta["neurons"] = ct.sum(axis=0)
meta.to_csv("results/animal_metadata.csv")

rds_genes = open("data/rds_genes.txt").read().split()
mats = []
for tag in ["Pman", "Ppol"]:
    a = sc.read_h5ad(f"data/{tag}_raw.h5ad")
    a = a[a.obs_names.isin(obs.index)].copy()
    nm = a.var["gene_name"].astype(str)
    keep = nm.isin(set(rds_genes)).values & ~nm.duplicated(keep="first").values
    a = a[:, keep].copy()
    a.var_names = nm[keep].values
    a = a[:, rds_genes].copy()
    mats.append((a.obs_names, a.X))
ids = np.concatenate([m[0] for m in mats])
X = sp.vstack([m[1] for m in mats]).tocsr()
obs = obs.loc[ids]
print("neurons exported:", X.shape)

m = pd.read_csv("data/matched_nuclei.csv").set_index("index")


def pseudobulk(labels, tag):
    keep = labels.notna().values
    groups = (labels[keep].astype(str) + "|" + obs.loc[keep, "animal"]).values
    codes, uniq = pd.factorize(groups)
    G = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
    PB = (X[keep].T @ G).tocsc()
    scipy.io.mmwrite(f"data/pseudobulk_counts{tag}.mtx", PB)
    pb_meta = pd.DataFrame({"group": uniq})
    pb_meta[["cluster", "animal"]] = pb_meta["group"].str.split("|", expand=True)
    pb_meta = pb_meta.join(meta[["species", "sex", "rep"]], on="animal")
    pb_meta["n_nuclei"] = np.asarray(G.sum(axis=0)).ravel()
    pb_meta.to_csv(f"data/pseudobulk_meta{tag}.csv", index=False)
    print("pseudobulk", tag or "(my clusters)", PB.shape)


pseudobulk(obs["cluster"], "")
pseudobulk(m["cluster_rds"].reindex(obs.index), "_auth")
pd.Series(rds_genes).to_csv("data/pseudobulk_genes.csv", index=False, header=False)


def write_pb(keep_mask, tag):
    labels = obs["cluster"].astype(str).values[keep_mask] + "|" + obs["animal"].values[keep_mask]
    codes, uniq = pd.factorize(labels)
    G = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
    PB = (X[keep_mask].T @ G).tocsc()
    scipy.io.mmwrite(f"data/pseudobulk_counts{tag}.mtx", PB)
    pm = pd.DataFrame({"group": uniq})
    pm[["cluster", "animal"]] = pm["group"].str.split("|", expand=True)
    pm = pm.join(meta[["species", "sex", "rep"]], on="animal")
    pm["n_nuclei"] = np.asarray(G.sum(axis=0)).ravel()
    pm.to_csv(f"data/pseudobulk_meta{tag}.csv", index=False)


grp = obs.groupby(["cluster", "animal"]).indices
cl_sp = obs.groupby(["cluster", "species"]).size().unstack(fill_value=0)


def draw(rng, frac_of):
    keep = np.zeros(len(obs), bool)
    for (cl, an), pos in grp.items():
        sp_ = obs["species"].values[pos[0]]
        n = int(round(frac_of(cl, sp_) * len(pos)))
        keep[rng.choice(pos, n, replace=False)] = True
    return keep


def eq(cl, sp_):
    n = cl_sp.loc[cl]; return min(1.0, n.min() / n[sp_]) if n[sp_] > 0 else 0.0


for k in range(1, 6):
    write_pb(draw(np.random.default_rng(100 + k), eq), f"_sub_eq_{k}")
for k in range(1, 4):
    write_pb(draw(np.random.default_rng(200 + k), lambda c, s: 0.5), f"_sub_f50_{k}")
    write_pb(draw(np.random.default_rng(300 + k), lambda c, s: 0.25), f"_sub_f25_{k}")
print("written")
