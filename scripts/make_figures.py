import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcdefaults()
MAN, POL, GREY = "#D98A1F", "#2F4FA8", "#B8B8B8"
SEX_M, SEX_F = "#2A7F62", "#B5446E"
INK = "#222222"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#888888", "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})

da = pd.read_csv("results/differential_abundance_edger.csv")
NAMES = {"my_clusters": {8: "Avp/Oxt", 6: "Gal/Moxd1"},
         "authors_clusters": {34: "Avp", 40: "Oxt", 20: "Gal/Th", 9: "Gal/Moxd1"}}
fig, axes = plt.subplots(2, 2, figsize=(14, 7.5), sharey="col")
for r, (an, title) in enumerate([("my_clusters", "My clusters"), ("authors_clusters", "Authors' clusters (RDS labels)")]):
    for c, (term, lab) in enumerate([("speciesPpol", "log2 FC, P. polionotus vs P. maniculatus"),
                                     ("sexM", "log2 FC, male vs female")]):
        ax = axes[r, c]
        d = da[(da.analysis == an) & (da.term == term)].copy()
        d["cl"] = d["cluster"].astype(int); d = d.sort_values("cl")
        if term == "speciesPpol":
            col = np.where(d.FDR < .05, np.where(d.logFC > 0, POL, MAN), GREY)
        else:
            col = np.where(d.FDR < .05, np.where(d.logFC > 0, SEX_M, SEX_F), GREY)
        ax.bar(np.arange(len(d)), d.logFC, color=col, width=0.7)
        ax.axhline(0, color="#888888", lw=0.8)
        ax.set_xticks(np.arange(len(d))); ax.set_xticklabels(d.cl, fontsize=6, rotation=90)
        ax.set_ylabel(lab); ax.set_title(f"{title}: {'species' if c == 0 else 'sex'} effect", loc="left", fontsize=10)
        for k, name in NAMES[an].items():
            i = list(d.cl).index(k); y = d.logFC.iloc[i]
            ax.annotate(name, (i, y), xytext=(0, 6 if y >= 0 else -12), textcoords="offset points",
                        ha="center", fontsize=7, color=INK)
axes[1, 0].set_xlabel("Neuron cluster"); axes[1, 1].set_xlabel("Neuron cluster")
fig.text(0.5, 0.005, "Coloured bars: FDR (BH) < 0.05; grey: not significant. Model: ~ rep + species + sex, TMM, QL F test.",
         ha="center", fontsize=8, color="#555555")
plt.tight_layout(rect=(0, 0.02, 1, 1)); plt.savefig("figures/differential_abundance.png", dpi=150); plt.close()

rows = []
for f, mode in [("results/sex_de_summary.csv", "strict filters"),
                ("results/sex_de_summary_authorlike.csv", "authors-like filters")]:
    if os.path.exists(f):
        s = pd.read_csv(f); s["mode"] = mode; rows.append(s)
S = pd.concat(rows)
S["setting"] = S["analysis"].map({"my_clusters": "my clusters", "authors_clusters": "authors' clusters"}) + "\n" + S["mode"]
settings = list(dict.fromkeys(S["setting"]))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for a, (colname, ttl) in zip(ax, [("genes_male_biased", "Male-biased genes"), ("genes_female_biased", "Female-biased genes")]):
    x = np.arange(len(settings)); w = 0.36
    for j, (sp_, colr, nm) in enumerate([("P.man", MAN, "P. maniculatus"), ("P.pol", POL, "P. polionotus")]):
        v = [S[(S.setting == s) & (S.species == sp_)][colname].iloc[0] for s in settings]
        b = a.bar(x + (j - .5) * w, v, w, color=colr, label=nm)
        for xi, vi in zip(x + (j - .5) * w, v):
            a.text(xi, vi + 3, str(vi), ha="center", fontsize=8, color=INK)
    if colname == "genes_male_biased":
        for j, (v, colr) in enumerate([(194, MAN), (70, POL)]):
            a.hlines(v, -0.5, len(settings) - 0.5, colors=colr, linestyles=":", lw=1.2)
        a.text(len(settings) - 0.5, 196, "paper: 194", ha="right", va="bottom", fontsize=8, color=MAN)
        a.text(len(settings) - 0.5, 72, "paper: 70", ha="right", va="bottom", fontsize=8, color=POL)
    a.set_xticks(x); a.set_xticklabels(settings, fontsize=8); a.set_title(ttl, loc="left", fontsize=10)
ax[0].set_ylabel("Genes with FDR < 0.05\n(counted once, at the cluster of highest logCPM)")
ax[1].legend(frameon=False, loc="upper right", fontsize=8)
plt.tight_layout(); plt.savefig("figures/sex_biased_genes.png", dpi=150); plt.close()
print(S.to_string(index=False))

plt.rcdefaults()
INK = "#222222"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888888", "axes.labelcolor": INK})
um = pd.read_csv("data/umap_neurons.csv", index_col=0)
lab = pd.read_csv("data/moffitt_labels_per_nucleus.csv", index_col=0)
tab = pd.read_csv("results/moffitt_cluster_labels_15pct.csv")
frac = pd.read_csv("results/moffitt_cluster_fractions.csv", index_col=0)

d = um.join(lab[["cluster", "neuron_type"]])
col = {"Excitatory": "#C8553D", "Inhibitory": "#3B6FB6"}
fig, ax = plt.subplots(figsize=(9, 8))
ax.scatter(d.UMAP_1, d.UMAP_2, s=0.4, c=d.neuron_type.map(col), alpha=0.5, linewidths=0, rasterized=True)
for cl, g in d.groupby("cluster"):
    r = tab[tab.cluster == cl].iloc[0]
    txt = r.top_label if r.top_fraction >= 0.15 else f"({r.top_label})"
    ax.text(g.UMAP_1.median(), g.UMAP_2.median(), f"{cl}: {txt}", fontsize=6.5, ha="center", va="center", color=INK,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
for k, v in col.items():
    ax.scatter([], [], c=v, s=25, label=k)
ax.legend(frameon=False, loc="lower left"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_title("Neuron clusters annotated with the best-matching Moffitt 2018 cluster\n(label shown as: my cluster: Moffitt cluster; brackets = below 15% rule)", fontsize=9, loc="left")
plt.tight_layout(); plt.savefig("figures/umap_moffitt_labels.png", dpi=150); plt.close()

cols = sorted(frac.columns, key=lambda s: (s[0] != "I", int(s.split("-")[1])))
order = tab.assign(k=tab.top_label.map({c: i for i, c in enumerate(cols)})).sort_values(["k", "top_fraction"], ascending=[True, False]).cluster
h = frac.loc[order.values, cols]
fig, ax = plt.subplots(figsize=(13, 9))
im = ax.imshow(h.values, aspect="auto", cmap="Blues", vmin=0, vmax=0.6)
ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=7)
ax.set_yticks(range(len(h))); ax.set_yticklabels([f"cluster {c}" for c in h.index], fontsize=7)
ax.axvline(sum(c.startswith("I") for c in cols) - 0.5, color=INK, lw=0.8)
ax.set_xlabel("Moffitt et al. 2018 neuron cluster (inhibitory I-, excitatory E-)"); ax.set_ylabel("My neuron cluster")
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01); cb.set_label("fraction of nuclei in my cluster (15% rule threshold = 0.15)")
cb.ax.axhline(0.15, color="#C8553D", lw=1)
plt.tight_layout(); plt.savefig("figures/cluster_vs_moffitt.png", dpi=150); plt.close()

m = pd.read_csv("data/matched_nuclei.csv").set_index("index")
x = lab.join(m["cluster_rds"], how="inner")
a = pd.crosstab(x["cluster_rds"], x["moffitt"], normalize="index")
out = pd.DataFrame({"top_label": a.idxmax(axis=1), "top_fraction": a.max(axis=1).round(2)})
out.index.name = "authors_cluster"; out.to_csv("results/authors_cluster_moffitt_top.csv")
print(out.loc[[20, 34, 40, 9, 5, 31]])

plt.rcdefaults()
rows = []
full = pd.read_csv("results/sex_de_summary.csv"); full = full[full.analysis == "my_clusters"]
for _, r in full.iterrows():
    rows.append(dict(cond="full data", draw=0, method="edgeR", species=r.species, male=r.genes_male_biased, female=r.genes_female_biased))
for f in sorted(glob.glob("results/robustness/edger_summary_sub_*.csv")):
    cond, k = re.search(r"sub_(eq|f50|f25)_(\d)", f).groups()
    name = {"eq": "equal nuclei\nper cluster", "f50": "50% of nuclei", "f25": "25% of nuclei"}[cond]
    for _, r in pd.read_csv(f).iterrows():
        rows.append(dict(cond=name, draw=int(k), method="edgeR", species=r.species, male=r.genes_male_biased, female=r.genes_female_biased))
for _, r in pd.read_csv("results/robustness/pydeseq2_summary.csv").iterrows():
    rows.append(dict(cond="full data", draw=0, method="PyDESeq2", species=r.species, male=r.genes_male_biased, female=r.genes_female_biased))
d = pd.DataFrame(rows)
w = d.pivot_table(index=["method", "cond", "draw"], columns="species", values="male").reset_index()
w["ratio_man_over_pol"] = (w["P.man"] / w["P.pol"]).round(2)
w.to_csv("results/robustness/robustness_table.csv", index=False)
print(w.to_string(index=False))
print(w[w.method == "edgeR"].groupby("cond")[["P.man", "P.pol", "ratio_man_over_pol"]].agg(["mean", "min", "max"]).round(2).to_string())

MAN, POL, INK = "#D98A1F", "#2F4FA8", "#222222"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
conds = ["full data", "equal nuclei\nper cluster", "50% of nuclei", "25% of nuclei"]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
e = w[w.method == "edgeR"]
x = np.arange(len(conds))
for j, (sp_, colr, nm) in enumerate([("P.man", MAN, "P. maniculatus"), ("P.pol", POL, "P. polionotus")]):
    for i, c in enumerate(conds):
        v = e[e.cond == c][sp_].values
        ax[0].scatter(np.full(len(v), i + (j - .5) * .25), v, color=colr, s=22, label=nm if i == 0 else None, zorder=3)
        ax[0].hlines(v.mean(), i + (j - .5) * .25 - .1, i + (j - .5) * .25 + .1, color=INK, lw=1.2, zorder=4)
ax[0].set_xticks(x); ax[0].set_xticklabels(conds, fontsize=8); ax[0].set_ylim(0, None)
ax[0].set_ylabel("Male-biased genes (edgeR, FDR < 0.05)"); ax[0].set_title("A. Counts under subsampling", loc="left", fontsize=10)
ax[0].legend(frameon=False, fontsize=8)
for i, c in enumerate(conds):
    v = e[e.cond == c]["ratio_man_over_pol"].values
    ax[1].scatter(np.full(len(v), i), v, color="#444444", s=22, zorder=3)
    ax[1].hlines(v.mean(), i - .15, i + .15, color=INK, lw=1.2, zorder=4)
ax[1].axhline(194 / 70, color="#888888", ls=":", lw=1.2); ax[1].text(3.4, 194 / 70 + .05, "paper 2.8", ha="right", fontsize=8, color="#555555")
ax[1].axhline(1, color="#BBBBBB", lw=0.8)
ax[1].set_xticks(x); ax[1].set_xticklabels(conds, fontsize=8); ax[1].set_ylim(0, None)
ax[1].set_ylabel("Ratio male-biased genes, P. man / P. pol"); ax[1].set_title("B. Ratio under subsampling", loc="left", fontsize=10)
m = w[(w.cond == "full data")].set_index("method")
labs = ["edgeR", "PyDESeq2"]; xx = np.arange(2); ww = .36
for j, (sp_, colr) in enumerate([("P.man", MAN), ("P.pol", POL)]):
    vals = [m.loc[l, sp_] for l in labs]
    b = ax[2].bar(xx + (j - .5) * ww, vals, ww, color=colr)
    for xi, v in zip(xx + (j - .5) * ww, vals): ax[2].text(xi, v + 1, int(v), ha="center", fontsize=8, color=INK)
ax[2].set_xticks(xx); ax[2].set_xticklabels([f"{l}\nratio {m.loc[l, 'ratio_man_over_pol']}" for l in labs])
ax[2].set_ylabel("Male-biased genes"); ax[2].set_title("C. Method: edgeR vs PyDESeq2 (full data)", loc="left", fontsize=10)
plt.tight_layout(); plt.savefig("figures/robustness_sex_bias.png", dpi=150)
