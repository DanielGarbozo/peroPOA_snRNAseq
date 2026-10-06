import sys
import numpy as np
import pandas as pd
import scipy.io
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_NUCLEI, FDR_CUT = 10, 0.05
np.random.seed(0)
tag = sys.argv[1] if len(sys.argv) > 1 else ""
out_name = sys.argv[2] if len(sys.argv) > 2 else "pydeseq2"
PB = scipy.io.mmread(f"data/pseudobulk_counts{tag}.mtx").tocsc()
genes = np.array(open("data/pseudobulk_genes.csv").read().split())
meta = pd.read_csv(f"data/pseudobulk_meta{tag}.csv")


def filter_by_expr(counts, sex, min_count=10, min_total=15):
    lib = counts.sum(axis=1)
    cpm_cut = min_count / np.median(lib) * 1e6
    n = min((sex == "F").sum(), (sex == "M").sum())
    cpm = counts.values / lib.values[:, None] * 1e6
    return ((cpm >= cpm_cut).sum(axis=0) >= n) & (counts.sum(axis=0).values >= min_total)


rows = []
for sp_ in ["P.man", "P.pol"]:
    for cl in sorted(meta["cluster"].unique(), key=int):
        m = meta[(meta.species == sp_) & (meta.cluster == cl) & (meta.n_nuclei >= MIN_NUCLEI)]
        if (m.sex == "F").sum() < 2 or (m.sex == "M").sum() < 2:
            continue
        n_par = 1 + (m.rep.nunique() - 1) + 1
        if len(m) - n_par < 2:
            continue
        cnt = pd.DataFrame(PB[:, m.index.values].T.toarray().astype(int), index=m["group"].values, columns=genes)
        keep = filter_by_expr(cnt, m["sex"].values)
        cnt = cnt.loc[:, keep]
        md = pd.DataFrame({"rep": m["rep"].astype(str).values, "sex": m["sex"].values}, index=m["group"].values)
        try:
            dds = DeseqDataSet(counts=cnt, metadata=md, design="~ rep + sex", quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["sex", "M", "F"], quiet=True)
            st.summary()
        except Exception as e:
            print(sp_, cl, "failed:", str(e)[:80], flush=True); continue
        r = st.results_df.dropna(subset=["padj"]).reset_index(names="gene")
        r["cluster"], r["species"] = cl, sp_
        rows.append(r)
    print(sp_, "done", flush=True)

res = pd.concat(rows)
res.to_csv(f"data/{out_name}_all.csv.gz", index=False)
sig = res[res.padj < FDR_CUT]
sig.to_csv(f"results/robustness/{out_name}_sig.csv", index=False)
summ = []
for sp_, d in sig.groupby("species"):
    top = d.loc[d.groupby("gene")["baseMean"].idxmax()]
    summ.append({"method": "PyDESeq2", "set": out_name, "species": sp_,
                 "clusters_tested": res[res.species == sp_].cluster.nunique(),
                 "sig_gene_cluster_pairs": len(d), "genes_male_biased": int((top.log2FoldChange > 0).sum()),
                 "genes_female_biased": int((top.log2FoldChange < 0).sum())})
summ = pd.DataFrame(summ); summ.to_csv(f"results/robustness/{out_name}_summary.csv", index=False)
print(summ.to_string(index=False))
