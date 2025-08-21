#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import sys
import tempfile
import re
from typing import List

# Headless plotting for HPC
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np


def atomic_write_csv(df: pd.DataFrame, out_path: Path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(out_path.parent), suffix=".csv") as tmp:
        tmp_name = tmp.name
        df.to_csv(tmp_name, index=False)
    os.replace(tmp_name, out_path)


def atomic_savefig(fig, out_path: Path, dpi=300):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", delete=False, dir=str(out_path.parent), suffix=".png") as tmp:
        tmp_name = tmp.name
        fig.savefig(tmp_name, dpi=dpi, bbox_inches="tight")
    os.replace(tmp_name, out_path)


def canon_id(s: str) -> str:
    # Strip literal [] and anything after the first space
    return re.sub(r'\[\]', '', str(s)).split()[0]


def compute_merged_coverage(intervals: List[tuple[int, int]], query_length: int) -> float:
    if not intervals:
        return 0.0
    # Sort intervals by start
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for current in intervals:
        if not merged or merged[-1][1] < current[0]:
            merged.append(current)
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], current[1]))
    union_length = sum(end - start + 1 for start, end in merged)
    return (union_length / query_length) * 100 if query_length > 0 else 0.0


def filter_invalid_lengths(blast_df: pd.DataFrame, min_length: int = 48, invalid_threshold: float = 0.5) -> pd.DataFrame:
    """Filter out rows with suspiciously short or missing query_length, log warnings."""
    valid_mask = blast_df['query_length'].fillna(0) >= min_length
    invalid_count = (~valid_mask).sum()
    if invalid_count > 0:
        print(f"WARNING: {invalid_count} rows with suspiciously short or missing query_length (<{min_length} bp) filtered out.")
        if invalid_count / len(blast_df) > invalid_threshold:
            raise ValueError(f"More than {invalid_threshold*100}% rows have invalid query_length; check input TSV and FASTA.")
    return blast_df[valid_mask].copy()


def calculate_coverage_stats(blast_file, fasta_file, output_dir):
    blast_file = Path(blast_file)
    fasta_file = Path(fasta_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_name = blast_file.stem  # e.g., ES01RB.unmapped_1_unclassified_fungi_90pct

    # Write files directly inside output_dir (no extra nesting)
    stats_csv = output_dir / f"{base_name}_coverage_stats.csv"
    det_csv   = output_dir / f"{base_name}_detailed_alignments.csv"
    plot_png  = output_dir / f"{base_name}_coverage_analysis.png"

    # Skip if all outputs exist and are non-empty
    if all(p.exists() and p.stat().st_size > 0 for p in (stats_csv, det_csv, plot_png)):
        print(f"Outputs already exist in '{output_dir}', skipping.")
        return None

    # Validate inputs
    if not blast_file.exists() or blast_file.stat().st_size == 0:
        print(f"ERROR: BLAST file missing or empty: {blast_file}", file=sys.stderr)
        sys.exit(2)
    if not fasta_file.exists() or fasta_file.stat().st_size == 0:
        print(f"ERROR: FASTA file missing or empty: {fasta_file}", file=sys.stderr)
        sys.exit(2)

    # Read BLAST results (updated: now with 13 columns, including appended query_length from phase 1)
    cols = [
        'qseqid','sseqid','pident','length','mismatch','gapopen',
        'qstart','qend','sstart','send','evalue','bitscore','query_length'
    ]
    blast_df = pd.read_csv(blast_file, sep='\t', header=None, names=cols)

    # Canonical FASTA IDs for unaligned contigs (still need FASTA for complete list)
    seq_lengths = {canon_id(r.id): len(r.seq) for r in SeqIO.parse(str(fasta_file), 'fasta')}

    # Debug: Print sample seq_lengths from FASTA
    print(f"Sample seq_lengths from FASTA: {dict(list(seq_lengths.items())[:5])}")  # First 5 for log

    if blast_df.empty:
        # Produce consistent empty outputs but still list FASTA contigs at 0% coverage
        print(f"WARNING: No rows in BLAST file: {blast_file}")
        all_fasta_ids = list(seq_lengths.keys())
        empty_stats = pd.DataFrame({
            "qseqid": all_fasta_ids,
            "pident_mean": 0, "pident_max": 0, "pident_std": 0,
            "query_coverage_mean": 0, "query_coverage_max": 0, "query_coverage_std": 0,
            "query_coverage_merged": 0,
            "evalue_min": np.nan, "bitscore_max": 0,
            "query_length": [seq_lengths[i] for i in all_fasta_ids]
        })
        atomic_write_csv(empty_stats, stats_csv)
        atomic_write_csv(blast_df, det_csv)
        fig, ax = plt.subplots(figsize=(6,4))
        ax.text(0.5, 0.5, "No BLAST alignments", ha="center", va="center")
        ax.axis("off")
        atomic_savefig(fig, plot_png)
        plt.close(fig)
        return empty_stats

    # Canonicalize qseqid (no need to remap query_length—it's already in the input)
    blast_df['cleaned_qseqid'] = blast_df['qseqid'].apply(canon_id)

    # Filter rows with invalid lengths (using the appended query_length column)
    blast_df = filter_invalid_lengths(blast_df)

    # If everything filtered out, still write consistent outputs with 0%-coverage for all FASTA contigs
    if blast_df.empty:
        print(f"WARNING: All BLAST rows filtered after length validation: {blast_file}")
        all_fasta_ids = list(seq_lengths.keys())
        empty_stats = pd.DataFrame({
            "qseqid": all_fasta_ids,
            "pident_mean": 0, "pident_max": 0, "pident_std": 0,
            "query_coverage_mean": 0, "query_coverage_max": 0, "query_coverage_std": 0,
            "query_coverage_merged": 0,
            "evalue_min": np.nan, "bitscore_max": 0,
            "query_length": [seq_lengths[i] for i in all_fasta_ids]
        })
        atomic_write_csv(empty_stats, stats_csv)
        atomic_write_csv(pd.DataFrame(columns=cols), det_csv)  # Updated cols
        fig, ax = plt.subplots(figsize=(6,4))
        ax.text(0.5, 0.5, "No BLAST alignments after filtering", ha="center", va="center")
        ax.axis("off")
        atomic_savefig(fig, plot_png)
        plt.close(fig)
        return empty_stats

    # Compute per-hit query coverage (%)
    blast_df['query_coverage'] = (blast_df['qend'] - blast_df['qstart'] + 1) / blast_df['query_length'] * 100

    # Aggregate per original qseqid for outputs (keep original qseqid label)
    stats = blast_df.groupby('qseqid').agg({
        'pident': ['mean','max','std'],
        'query_coverage': ['mean','max','std'],
        'evalue': 'min',
        'bitscore': 'max'
    }).round(2)
    stats.columns = ['_'.join(c) for c in stats.columns]
    stats = stats.reset_index()

    # Provide a unique query_length per qseqid (e.g., max; lengths per query should be identical)
    unique_query_lengths = blast_df.groupby('qseqid')['query_length'].max()
    stats['query_length'] = stats['qseqid'].map(unique_query_lengths)

    # Compute merged coverage per query
    grouped = blast_df.groupby('qseqid')
    merged_coverages = {}
    for qid, group in grouped:
        intervals = list(zip(group['qstart'], group['qend']))
        qlen = group['query_length'].max()  # Assuming consistent length
        merged_coverages[qid] = compute_merged_coverage(intervals, qlen)
    stats['query_coverage_merged'] = stats['qseqid'].map(merged_coverages)

    # Optionally add 0%-coverage rows for FASTA contigs with no alignments
    # (helps downstream by ensuring all contigs appear in stats)
    all_fasta_ids = set(seq_lengths.keys())
    aligned_ids = set(stats['qseqid'].apply(canon_id))
    unaligned_ids = all_fasta_ids - aligned_ids
    if unaligned_ids:
        unaligned_stats = pd.DataFrame({
            "qseqid": list(unaligned_ids),
            "pident_mean": 0, "pident_max": 0, "pident_std": 0,
            "query_coverage_mean": 0, "query_coverage_max": 0, "query_coverage_std": 0,
            "query_coverage_merged": 0,
            "evalue_min": np.nan, "bitscore_max": 0,
            "query_length": [seq_lengths[i] for i in unaligned_ids]
        })
        stats = pd.concat([stats, unaligned_stats], ignore_index=True)

    # Standardize all qseqid to canonical form for consistency
    stats['qseqid'] = stats['qseqid'].apply(canon_id)

    # Save outputs
    atomic_write_csv(stats, stats_csv)
    atomic_write_csv(blast_df, det_csv)

    # Plot (updated to use merged coverage)
    plot_coverage_analysis(blast_df, stats, plot_png)
    return stats


def plot_coverage_analysis(blast_df, stats, out_png_path):
    fig, ax = plt.subplots(2, 2, figsize=(12, 10))

    # 1) Query Coverage distribution (per-hit)
    ax[0,0].hist(blast_df['query_coverage'], bins=50, color='skyblue', edgecolor='black')
    ax[0,0].set(title='Per-Hit Query Coverage Distribution', xlabel='Coverage (%)', ylabel='Count')

    # 2) Identity vs Coverage (per-hit)
    ax[0,1].scatter(blast_df['query_coverage'], blast_df['pident'], alpha=0.5)
    ax[0,1].set(title='Identity vs Per-Hit Coverage', xlabel='Coverage (%)', ylabel='Percent Identity')

    # 3) Query length distribution
    ax[1,0].hist(stats['query_length'], bins=50, color='lightgreen', edgecolor='black')
    ax[1,0].set(title='Query Length Distribution', xlabel='Length (bp)', ylabel='Count')

    # 4) Merged coverage vs query length
    ax[1,1].scatter(stats['query_length'], stats['query_coverage_merged'], alpha=0.5)
    ax[1,1].set(title='Merged Coverage vs Query Length', xlabel='Length (bp)', ylabel='Merged Coverage (%)')

    plt.tight_layout()
    atomic_savefig(fig, Path(out_png_path))
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='BLAST coverage analysis')
    p.add_argument('--blast',      required=True, help='BLAST TSV file')
    p.add_argument('--fasta',      required=True, help='Query FASTA file')
    p.add_argument('--output_dir', required=True, help='Output directory (files will be saved with auto-generated names)')
    args = p.parse_args()
    calculate_coverage_stats(args.blast, args.fasta, args.output_dir)
