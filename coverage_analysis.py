#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import re
from pathlib import Path
from typing import List, Tuple, Dict
from Bio import SeqIO

def canon_id(s: str) -> str:
    return re.sub(r'\[\]', '', str(s)).split()[0]

def calculate_merged_coverage(intervals: List[Tuple[int, int]], query_length: int) -> float:
    """
    Calculate merged % coverage (unique bases covered / query_length * 100).
    """
    if query_length <= 0 or not intervals:
        return 0.0
    coverage_array = np.zeros(query_length, dtype=bool)
    for start, end in intervals:
        start_idx = max(0, start - 1)
        end_idx = min(query_length, end)
        coverage_array[start_idx:end_idx] = True
    covered_bases = np.count_nonzero(coverage_array)
    return (covered_bases / query_length) * 100

def process_tsv(tsv_path: Path, fasta_dir: Path, output_dir: Path, coverage_threshold: float = 70.0):
    # Read TSV
    blast_df = pd.read_csv(tsv_path, sep='\t', names=[
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ])

    # Extract sample and db
    sample = tsv_path.stem.split('_vs_')[0]
    db = tsv_path.stem.split('_vs_')[1].replace('_filtered', '')

    # Load FASTA
    fasta_path = fasta_dir / f"{sample}.clean.fasta"
    if not fasta_path.exists():
        print(f"SKIP: Missing FASTA for {sample}: {fasta_path}")
        return

    seq_records = list(SeqIO.parse(str(fasta_path), 'fasta'))
    seq_lengths = {canon_id(r.id): len(r.seq) for r in seq_records}
    print(f"Loaded {len(seq_lengths)} sequences for {sample} from {fasta_path}.")

    if blast_df.empty:
        print(f"WARNING: Empty TSV: {tsv_path}")
        results = pd.DataFrame({
            'qseqid': list(seq_lengths.keys()),
            'num_alignments': 0,
            'pident_mean': 0.0,
            'coverage_pct': 0.0,
            'query_length': list(seq_lengths.values())
        })
    else:
        blast_df['cleaned_qseqid'] = blast_df['qseqid'].apply(canon_id)
        blast_df = blast_df[blast_df['cleaned_qseqid'].isin(seq_lengths)]

        grouped = blast_df.groupby('cleaned_qseqid')
        results_list = []
        for qid, group in grouped:
            intervals = list(zip(group['qstart'], group['qend']))
            qlen = seq_lengths[qid]
            coverage = calculate_merged_coverage(intervals, qlen)
            results_list.append({
                'qseqid': qid,
                'num_alignments': len(group),
                'pident_mean': group['pident'].mean(),
                'coverage_pct': coverage,
                'query_length': qlen
            })

        # Add unaligned queries
        aligned_ids = {r['qseqid'] for r in results_list}
        for qid, qlen in seq_lengths.items():
            if qid not in aligned_ids:
                results_list.append({
                    'qseqid': qid,
                    'num_alignments': 0,
                    'pident_mean': 0.0,
                    'coverage_pct': 0.0,
                    'query_length': qlen
                })

        results = pd.DataFrame(results_list)

    # Save full results
    out_csv = output_dir / f"{sample}_vs_{db}_coverage_stats.csv"
    results.to_csv(out_csv, index=False)
    print(f"Saved full stats to {out_csv}")

    # Filter and extract high-coverage queries
    high_cov = results[results['coverage_pct'] > coverage_threshold]
    if not high_cov.empty:
        high_cov.to_csv(output_dir / f"{sample}_vs_{db}_high_coverage.csv", index=False)
        print(f"Extracted {len(high_cov)} queries with >{coverage_threshold}% coverage to high_coverage CSV.")

        # Extract sequences to new FASTA
        high_ids = set(high_cov['qseqid'])
        high_fasta = output_dir / f"{sample}_vs_{db}_high_coverage.fasta"
        with open(high_fasta, 'w') as f:
            for record in seq_records:
                if canon_id(record.id) in high_ids:
                    SeqIO.write(record, f, 'fasta')
        print(f"Extracted high-coverage sequences to {high_fasta}")
    else:
        print(f"No queries with >{coverage_threshold}% coverage.")

def main(tsv_dir: str, fasta_dir: str, output_dir: str, coverage_threshold: float):
    tsv_dir = Path(tsv_dir)
    fasta_dir = Path(fasta_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for tsv_path in sorted(tsv_dir.glob('*_filtered.tsv')):
        process_tsv(tsv_path, fasta_dir, output_dir, coverage_threshold)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate coverage % and extract >70% queries from filtered BLAST TSVs.')
    parser.add_argument('--tsv_dir', default='../filtered_blast_results', help='Directory with filtered TSVs')
    parser.add_argument('--fasta_dir', default='../cleaned_fastas', help='Directory with clean FASTAs')
    parser.add_argument('--output_dir', default='coverage_results', help='Output directory')
    parser.add_argument('--coverage_threshold', type=float, default=70.0, help='Coverage threshold for extraction')
    args = parser.parse_args()
    main(args.tsv_dir, args.fasta_dir, args.output_dir, args.coverage_threshold)
