import gzip
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import scrublet as scr
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.random.seed(0)
RAW = "data/raw/"
GSM = {"rep01": "GSM8409692", "rep02": "GSM8409693", "rep03": "GSM8409694",
       "rep04": "GSM8409695", "rep05": "GSM8409696", "rep06": "GSM8409697"}
SPECIES = {"P.man": "p_man", "P.pol": "p_pol"}
NAMEMAP = {"P.man": "data/GSE272719_Pman.gene_name.txt.gz",
           "P.pol": "data/GSE272719_Ppol.gene_name.txt.gz"}


def read_mtx_columns(path, keep_cols, n_cols_total):
    chunks = pd.read_csv(path, sep=" ", skiprows=3, header=None,
                         names=["gene", "cell", "count"], chunksize=5_000_000,
                         dtype={"gene": np.int32, "cell": np.int32, "count": np.int32})
    parts = [c[c["cell"].isin(keep_cols)] for c in chunks]
    return pd.concat(parts)


adatas = {sp_: [] for sp_ in SPECIES}
for rep, gsm in GSM.items():
    demux = pd.read_csv(f"{RAW}{gsm}_demultiplexed_sample_barcodes.{rep}.txt.gz", sep="\t")
    for species, tag in SPECIES.items():
        sub = demux[demux["species"] == species].copy()
        barcodes = pd.read_csv(f"{RAW}{gsm}_barcodes.{rep}.{tag}.tsv.gz", header=None)[0]
        feats = pd.read_csv(f"{RAW}{gsm}_features.{rep}.{tag}.tsv.gz", sep="\t", header=None)
        with gzip.open(f"{RAW}{gsm}_matrix.{rep}.{tag}.mtx.gz", "rt") as fh:
            fh.readline(); fh.readline()
            n_genes, n_cells, _ = map(int, fh.readline().split())
        pos = pd.Series(np.arange(1, len(barcodes) + 1), index=barcodes)
        sub = sub[sub["barcode"].isin(pos.index)]
        sub["col"] = pos.loc[sub["barcode"]].values
        trip = read_mtx_columns(f"{RAW}{gsm}_matrix.{rep}.{tag}.mtx.gz", set(sub["col"]), n_cells)
        remap = pd.Series(np.arange(len(sub)), index=sub["col"].values)
        X = sp.csr_matrix((trip["count"].values,
                           (remap.loc[trip["cell"]].values, trip["gene"].values - 1)),
                          shape=(len(sub), n_genes))
        a = ad.AnnData(X)
        a.obs_names = (rep + "_" + sub["barcode"]).values
        a.obs["rep"] = rep
        a.obs["species"] = species
        a.obs["sex"] = sub["sex"].str.upper().values
        a.obs["barcode"] = sub["barcode"].values
        a.var_names = feats[0].values
        adatas[species].append(a)
        print(rep, species, a.shape, flush=True)

for species, lst in adatas.items():
    a = ad.AnnData(sp.vstack([x.X for x in lst]).tocsr(),
                   obs=pd.concat([x.obs for x in lst]))
    a.var_names = lst[0].var_names
    nm = pd.read_csv(NAMEMAP[species], sep="\t", comment=None)
    nm.columns = ["gene_id", "gene_name"]
    nm = nm.drop_duplicates("gene_id").set_index("gene_id")["gene_name"]
    a.var["gene_name"] = nm.reindex(a.var_names).values
    a.obs["animal"] = a.obs["rep"] + "_" + a.obs["species"] + "_" + a.obs["sex"]
    a.var["mt"] = a.var["gene_name"].astype(str).str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    a.write_h5ad(f"data/{species.replace('.', '')}_raw.h5ad")
    print(species, a.shape, "mt genes:", int(a.var["mt"].sum()))

np.random.seed(0)
sc.settings.verbosity = 1
PAPER_TOTAL = 95057

tables, obs_all = [], []
for species in ["Pman", "Ppol"]:
    a = sc.read_h5ad(f"data/{species}_raw.h5ad")
    pat = re.compile(r"^(mt-|MT-|Nd[1-6]l?$|Cox[123]$|Atp[68]$|Cytb$)", re.I)
    mito = a.var_names[a.var["gene_name"].astype(str).str.match(pat)]
    print(species, "mitochondrial-like genes in annotation:", len(mito), list(mito[:5]))

    xist = a.var_names[a.var["gene_name"] == "Xist"]
    Xc = a[:, xist].X.toarray().ravel()
    a.obs["xist_counts"] = Xc
    a.obs["xist_cpm"] = Xc / a.obs["total_counts"].values * 1e6

    a.obs["doublet_score"] = np.nan
    a.obs["predicted_doublet"] = False
    for animal, idx in a.obs.groupby("animal").indices.items():
        sub = a.X[idx]
        s = scr.Scrublet(sub, expected_doublet_rate=0.06, random_state=0)
        score, pred = s.scrub_doublets(verbose=False, n_prin_comps=20)
        a.obs.iloc[idx, a.obs.columns.get_loc("doublet_score")] = score
        a.obs.iloc[idx, a.obs.columns.get_loc("predicted_doublet")] = pred
    a.obs["species_tag"] = species
    obs_all.append(a.obs.copy())

obs = pd.concat(obs_all)
obs.to_csv("data/qc_per_nucleus.csv")

tab = (obs.groupby(["rep", "species", "sex"]).agg(
    nuclei=("total_counts", "size"),
    median_umis=("total_counts", "median"),
    median_genes=("n_genes_by_counts", "median"),
    doublets=("predicted_doublet", "sum")).reset_index())
tab.to_csv("results/qc_nuclei_per_sample.csv", index=False)
print(tab.to_string(index=False))
print("TOTAL nuclei:", len(obs), "paper:", PAPER_TOTAL)

sex_chk = obs.groupby(["species", "sex"]).agg(
    frac_xist_pos=("xist_counts", lambda x: (x > 0).mean()),
    median_xist_cpm=("xist_cpm", "median")).reset_index()
sex_chk.to_csv("results/qc_sex_xist.csv", index=False)
print(sex_chk.to_string(index=False))

obs["sample"] = obs["rep"].astype(str) + "_" + obs["species"].astype(str).str.replace("P.", "", regex=False) + "_" + obs["sex"].astype(str)
fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
for a_, col, lab in zip(ax, ["total_counts", "n_genes_by_counts", "doublet_score"],
                        ["UMIs per nucleus (log10)", "Genes per nucleus (log10)", "Scrublet score"]):
    d = obs.assign(v=np.log10(obs[col]) if col != "doublet_score" else obs[col])
    sns.violinplot(data=d, x="sample", y="v", hue="species", ax=a_, cut=0, linewidth=.5, density_norm="width")
    a_.set_ylabel(lab); a_.legend(loc="upper right", fontsize=7)
plt.xticks(rotation=90, fontsize=7); plt.tight_layout()
plt.savefig("figures/supplemental/qc_per_sample.png", dpi=150)
