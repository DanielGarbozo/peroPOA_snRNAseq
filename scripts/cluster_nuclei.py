import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import harmonypy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.random.seed(0)
sc.settings.verbosity = 2
sc.settings.figdir = "figures/supplemental"

parts = {}
for s in ["Pman", "Ppol"]:
    a = sc.read_h5ad(f"data/{s}_raw.h5ad")
    a = a[:, a.var["gene_name"].notna()].copy()
    a = a[:, ~a.var["gene_name"].duplicated(keep="first")].copy()
    a.var_names = a.var["gene_name"].astype(str).values
    parts[s] = a
shared = sorted(set(parts["Pman"].var_names) & set(parts["Ppol"].var_names))
print("shared genes:", len(shared))
X = sp.vstack([parts[s][:, shared].X for s in ["Pman", "Ppol"]]).tocsr()
obs = pd.concat([parts[s].obs[["rep", "species", "sex", "animal", "barcode", "total_counts",
                               "n_genes_by_counts"]] for s in ["Pman", "Ppol"]])
adata = ad.AnnData(X, obs=obs)
adata.var_names = shared
adata.obs["species"] = adata.obs["species"].astype(str)
adata.obs["rep"] = adata.obs["rep"].astype(str)
q = pd.read_csv("data/qc_per_nucleus.csv", index_col=0)
adata.obs["doublet_score"] = q.loc[adata.obs_names, "doublet_score"].values
adata.obs["predicted_doublet"] = q.loc[adata.obs_names, "predicted_doublet"].values
sc.pp.filter_genes(adata, min_cells=3)
print(adata)
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="rep")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, random_state=0)

ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["rep", "species"], random_state=0)
Z = np.asarray(ho.Z_corr)
adata.obsm["X_pca_harmony"] = Z.T if Z.shape[0] == adata.obsm["X_pca"].shape[1] else Z
print("harmony shape:", adata.obsm["X_pca_harmony"].shape)

sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=30, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.tl.leiden(adata, resolution=0.8, random_state=0, flavor="igraph", n_iterations=2,
             directed=False, key_added="leiden_all")

MARKERS = {
    "Neuron": ["Snap25", "Syt1", "Rbfox3", "Stmn2", "Grin1"],
    "Astrocyte": ["Slc1a2", "Slc1a3", "Aqp4", "Gfap", "Aldh1l1"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mobp"],
    "OPC": ["Pdgfra", "Cspg4", "Vcan"],
    "Microglia": ["Csf1r", "Cx3cr1", "P2ry12", "Tmem119"],
    "Endothelial": ["Flt1", "Cldn5", "Pecam1"],
    "Mural/VLMC": ["Pdgfrb", "Acta2", "Myh11", "Ptgds", "Cfh", "Slc6a13", "Igfbp2", "Dcn"],
    "Ependymal": ["Foxj1", "Ccdc153"],
    "Macrophage": ["Mrc1", "F13a1"],
}
raw = adata.raw.to_adata()
for k, g in MARKERS.items():
    g = [x for x in g if x in raw.var_names]
    sc.tl.score_genes(raw, g, score_name="score_" + k, random_state=0)
sc_cols = ["score_" + k for k in MARKERS]
cl_score = raw.obs.groupby(adata.obs["leiden_all"], observed=True)[sc_cols].mean()
cl_score.columns = list(MARKERS)
glial = [c for c in cl_score.columns if c != "Neuron"]
glial_best = cl_score[glial].max(axis=1)
is_neuron = (cl_score["Neuron"] >= 0.5) & (cl_score["Neuron"] - glial_best >= 0.2)
cl_class = cl_score[glial].idxmax(axis=1).where(~is_neuron, "Neuron")
adata.obs["major_class"] = adata.obs["leiden_all"].map(cl_class).astype(str)
cl_score.assign(assigned=cl_class).to_csv("results/cluster_marker_scores.csv")
print(adata.obs["major_class"].value_counts())

adata.obs["is_neuron_cluster"] = adata.obs["major_class"] == "Neuron"

adata.write_h5ad("data/all_nuclei_integrated.h5ad")
raw.obs.to_csv("data/obs_all_nuclei.csv")

fig, ax = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.umap(adata, color="major_class", ax=ax[0], show=False, title="Major class", legend_fontsize=7)
sc.pl.umap(adata, color="species", ax=ax[1], show=False, title="Species")
sc.pl.umap(adata, color="rep", ax=ax[2], show=False, title="Replicate")
plt.tight_layout(); plt.savefig("figures/supplemental/umap_all_nuclei.png", dpi=130)
genes = [g for v in MARKERS.values() for g in v if g in raw.var_names]
sc.pl.dotplot(raw, genes, groupby="leiden_all", show=False, save="_markers.png")

np.random.seed(0)
TARGET = 53

full = sc.read_h5ad("data/all_nuclei_integrated.h5ad")
raw = full.raw.to_adata()
raw.obs = full.obs.copy()
neu = raw[raw.obs["is_neuron_cluster"].values].copy()
print("neurons:", neu.n_obs, "(paper: 52,121)")

sc.pp.highly_variable_genes(neu, n_top_genes=3000, flavor="seurat", batch_key="species")
neu.raw = neu
neu = neu[:, neu.var["highly_variable"]].copy()
sc.pp.scale(neu, max_value=10)
sc.tl.pca(neu, n_comps=50, random_state=0)
ho = harmonypy.run_harmony(neu.obsm["X_pca"], neu.obs, ["rep", "species"], theta=[2, 6],
                           max_iter_harmony=20, random_state=0)
Z = np.asarray(ho.Z_corr)
neu.obsm["X_pca_harmony"] = Z.T if Z.shape[0] == 50 else Z
sc.pp.neighbors(neu, use_rep="X_pca_harmony", n_neighbors=30, random_state=0)
sc.tl.umap(neu, random_state=0)

best = None
for res in [1.0, 2.2, 3.0, 3.5, 4.0, 5.0]:
    key = f"leiden_{res}"
    sc.tl.leiden(neu, resolution=res, random_state=0, flavor="igraph", n_iterations=2,
                 directed=False, key_added=key)
    n = neu.obs[key].nunique()
    print("resolution", res, "->", n, "clusters", flush=True)
    if best is None or abs(n - TARGET) < abs(best[1] - TARGET):
        best = (res, n)
print("chosen resolution:", best)
neu.obs["cluster"] = neu.obs[f"leiden_{best[0]}"].astype(str)

r = neu.raw.to_adata()
ex = np.asarray(r[:, "Slc17a6"].X.todense()).ravel()
inh = np.asarray(r[:, ["Gad1", "Gad2"]].X.todense()).mean(axis=1)
neu.obs["ex"], neu.obs["inh"] = ex, inh
cl = neu.obs.groupby("cluster")[["ex", "inh"]].mean()
cl["type"] = np.where(cl["ex"] > cl["inh"], "Excitatory", "Inhibitory")
neu.obs["neuron_type"] = neu.obs["cluster"].map(cl["type"])
cl.to_csv("results/neuron_cluster_ei.csv")
print(neu.obs["neuron_type"].value_counts())

neu.obs.to_csv("data/obs_neurons.csv")
pd.DataFrame(neu.obsm["X_umap"], index=neu.obs_names, columns=["UMAP_1", "UMAP_2"]).to_csv("data/umap_neurons.csv")
neu.write_h5ad("data/neurons_integrated.h5ad")

fig, ax = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.umap(neu, color="cluster", ax=ax[0], show=False, legend_loc="on data", legend_fontsize=6, title="Neuron clusters")
sc.pl.umap(neu, color="neuron_type", ax=ax[1], show=False, title="Excitatory / inhibitory")
sc.pl.umap(neu, color="species", ax=ax[2], show=False, title="Species")
plt.tight_layout(); plt.savefig("figures/supplemental/umap_neurons.png", dpi=130)
