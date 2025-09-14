#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Function to extract contig lengths from a single FASTA file
def extract_lengths_from_fasta(fasta_path):
    contigs = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Header line, e.g., >k127_3205 flag=1 multi=3.0000 len=583
                parts = line.strip().split()
                qseqid = parts[0][1:]  # Remove '>' to get k127_3205
                # Find len= part
                length = None
                for p in parts:
                    if p.startswith('len='):
                        length = int(p.split('=')[1])
                        break
                if length is not None:
                    contigs.append({'qseqid': qseqid, 'length_bp': length})
    return pd.DataFrame(contigs)

# Directory with cleaned FASTA files
fasta_dir = Path("../cleaned_fastas")

# List all .clean.fasta files
fasta_files = list(fasta_dir.glob("*.clean.fasta"))

# Extract lengths from all FASTA files, adding Sample column
lengths_dfs = []
for fasta_path in fasta_files:
    sample = fasta_path.stem.replace(".clean", "")  # e.g., HT02RT
    df = extract_lengths_from_fasta(fasta_path)
    df["Sample"] = sample
    lengths_dfs.append(df)

if not lengths_dfs:
    raise FileNotFoundError("No .clean.fasta files found in ../cleaned_fastas/")

lengths = pd.concat(lengths_dfs, ignore_index=True)

# Load coverage data
coverage_dir = Path("coverage_files")
coverage_files = list(coverage_dir.glob("*_capped_stats.csv"))
cov_dfs = []
for f in coverage_files:
    df = pd.read_csv(f)
    base = f.stem.replace("_capped_stats", "")
    parts = base.split("_vs_")[-1] if "_vs_" in base else base
    grp = "OTHER"
    for k, g in [("bact", "BACT"), ("virus", "VIRUS"), ("parasite", "PARASITE"), ("fungi", "FUNGI")]:
        if k in parts.lower():
            grp = g
            break
    df["MicrobialGroup"] = grp
    # Infer Sample from base, adjusting for possible mismatches (e.g., use prefix like HT02RB)
    sample = base.split("_vs_")[0] if "_vs_" in base else base
    df["Sample"] = sample
    cov_dfs.append(df)

if not cov_dfs:
    raise FileNotFoundError("No *_capped_stats.csv files found in coverage_files/")

cov = pd.concat(cov_dfs, ignore_index=True)

# Choose coverage column to use
possible_coverage_cols = ["capped_coverage_pct", "query_coverage_merged", "query_coverage_mean"]
coverage_col = None
for col in possible_coverage_cols:
    if col in cov.columns:
        coverage_col = col
        break
if coverage_col is None:
    raise ValueError(f"None of the expected coverage columns found: {possible_coverage_cols}")

# Merge coverage with lengths on Sample and qseqid
covL = cov.merge(lengths, on=["Sample", "qseqid"], how="left").dropna(subset=["length_bp", coverage_col])

# If merge results in empty DataFrame, print warning
if covL.empty:
    print("Warning: No matching Sample/qseqid pairs found between coverage data and FASTA lengths. Check sample naming consistency (e.g., HT02RB vs HT02RT).")

# Plot scatter with log-log scale
plt.figure(figsize=(10,7))
sns.scatterplot(data=covL, x="length_bp", y=coverage_col, hue="MicrobialGroup", size="num_alignments", sizes=(10, 200), palette="Set2", alpha=0.7, edgecolor=None)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Contig length (bp) [log scale]")
plt.ylabel(f"{coverage_col} [log scale]")
plt.title("Scatter plot of Contig length vs Coverage by Microbial Group (Multi-Sample)")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig("plots/contig_length_vs_coverage_scatter_multi.png", dpi=300)
plt.close()
