#!/usr/bin/env python3
"""
Comprehensive coverage analysis with robust scatter plots
- Guards against log-scale collapse (many values == 1)
- Adds jitter/auto linear switch for visibility
- Keeps informative sample-level bubble plot
- Adds alignment-count distribution plot
"""

#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Inputs (adjust paths to your compiled coverage stats)
coverage_csvs = list(Path("coverage_files").glob("coverage_*/*_coverage_stats.csv"))
dfs = []
for p in coverage_csvs:
    df = pd.read_csv(p)
    # Expected columns in your Phase2 stats: qseqid, query_coverage_mean, query_coverage_merged, bitscore_max, pident_mean, Sample (if you added), etc.
    # Infer Sample & MicrobialGroup from filename if not in CSV
    base = p.stem.replace("_coverage_stats","")
    # coverage_ES01RB.unmapped_1_classified_bact_90pct
    parts = base.split("coverage_")[-1]
    # crude group inference
    grp = "OTHER"
    for k,g in [("bact","BACT"),("virus","VIRUS"),("parasite","PARASITE"),("fungi","FUNGI")]:
        if k in parts.lower():
            grp=g; break
    df["MicrobialGroup"] = grp
    df["Sample"] = parts  # optional: preserve full name; later you can parse it down
    dfs.append(df)

cov = pd.concat(dfs, ignore_index=True)

# Choose coverage signal
y = "query_coverage_merged" if "query_coverage_merged" in cov.columns else "query_coverage_mean"
cov = cov[cov[y].notnull()].copy()

plt.figure(figsize=(12,7))
ax = sns.violinplot(data=cov, x="MicrobialGroup", y=y, inner=None, cut=0, linewidth=0.8,
                    order=["BACT","FUNGI","PARASITE","VIRUS"])
sns.boxplot(data=cov, x="MicrobialGroup", y=y, showcaps=True, boxprops={'facecolor':'white','alpha':0.7},
            showfliers=False, whiskerprops={'linewidth':1}, width=0.25,
            order=["BACT","FUNGI","PARASITE","VIRUS"], ax=ax)
ax.set_yscale("log")
ax.set_ylabel("Merged query coverage (%) [log]")
ax.set_xlabel("Microbial Group")
ax.set_title("Per-contig coverage distribution by microbial group")
ax.grid(True, axis='y', alpha=0.2)
plt.tight_layout()
Path("plots").mkdir(exist_ok=True)
plt.savefig("plots/violin_coverage_by_group.png", dpi=300)

import numpy as np

# cov from previous snippet
y = "query_coverage_merged" if "query_coverage_merged" in cov.columns else "query_coverage_mean"

# summarize counts and mean per sample×group
agg = (cov.groupby(["Sample","MicrobialGroup"], as_index=False)
         .agg(n_contigs=("qseqid","nunique"),
              mean_cov=(y,"mean")))

# Select a manageable set of samples (e.g., all ES* or top N by contigs)
# Example: keep all samples; for many samples, you can facet-wrap by group
g = sns.FacetGrid(agg, col="MicrobialGroup", col_wrap=2, sharey=False, height=4)
def _plot(data, color, **k):
    # box-ish: show distribution of per-sample means (not a real boxplot but conveys spread)
    sns.stripplot(data=data, x="Sample", y="mean_cov", color=color, alpha=0.6, size=5)
    # overlay bubbles sized by n_contigs
    s = np.clip(data["n_contigs"]*0.3, 10, 200)
    plt.scatter(x=np.arange(len(data["Sample"])), y=data["mean_cov"], s=s, facecolors="none", edgecolors="k", linewidths=0.6)
    plt.xticks(rotation=90)

g.map_dataframe(_plot)
g.set_axis_labels("Sample", "Mean coverage (%)")
for ax in g.axes.flat:
    ax.set_yscale("log")
    ax.grid(True, axis='y', alpha=0.2)

plt.tight_layout()
plt.savefig("plots/box_scatter_sample_group_coverage.png", dpi=300)

# contigs per sample×group
counts = (cov.groupby(["Sample","MicrobialGroup"], as_index=False)
            .agg(n_contigs=("qseqid","nunique")))

totals = counts.groupby("Sample", as_index=False)["n_contigs"].sum().rename(columns={"n_contigs":"total"})
counts = counts.merge(totals, on="Sample", how="left")
counts["rel"] = counts["n_contigs"] / counts["total"]

# order samples by total (or by category if you have it)
sample_order = totals.sort_values("total", ascending=False)["Sample"].tolist()
group_order = ["BACT","FUNGI","PARASITE","VIRUS"]

# pivot for stacked bars
pvt = counts.pivot(index="Sample", columns="MicrobialGroup", values="rel").fillna(0).reindex(sample_order)
pvt = pvt.reindex(columns=[g for g in group_order if g in pvt.columns])

colors = {"BACT":"#1f77b4","FUNGI":"#9467bd","PARASITE":"#2ca02c","VIRUS":"#ff7f0e"}
ax = pvt.plot(kind="bar", stacked=True, figsize=(14,6), color=[colors[c] for c in pvt.columns])
ax.set_ylabel("Relative abundance (by contig count)")
ax.set_xlabel("Sample")
ax.set_title("Relative abundance of microbial groups per sample")
ax.legend(title="MicrobialGroup", bbox_to_anchor=(1.02,1), loc="upper left")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("plots/stacked_relative_abundance_per_sample.png", dpi=300)

# Suppose you have a TSV mapping: contig_lengths.tsv with columns: qseqid, length_bp
lengths = pd.read_csv("contig_lengths.tsv", sep="\t")
covL = cov.merge(lengths, on="qseqid", how="left").dropna(subset=["length_bp"])
y = "query_coverage_merged" if "query_coverage_merged" in covL.columns else "query_coverage_mean"

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8,6))
sns.kdeplot(data=covL, x="length_bp", y=y, fill=True, thresh=0.05, levels=30, cmap="mako")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Contig length (bp) [log]")
plt.ylabel("Coverage (%) [log]")
plt.title("Density: Coverage vs contig length")
plt.tight_layout()
plt.savefig("plots/density_length_vs_coverage.png", dpi=300)
plt.close()

plt.figure(figsize=(8,6))
hb = plt.hexbin(covL["length_bp"], covL[y], gridsize=50, bins="log", cmap="viridis")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("Contig length (bp) [log]"); plt.ylabel("Coverage (%) [log]")
cb = plt.colorbar(hb); cb.set_label("log10(N)")
plt.title("Hexbin: Coverage vs contig length")
plt.tight_layout()
plt.savefig("plots/hexbin_length_vs_coverage.png", dpi=300)
plt.close()

