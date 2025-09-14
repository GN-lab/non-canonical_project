#!/bin/bash
#SBATCH --job-name=phase2
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8042
#SBATCH --time=24:00:00
#SBATCH --output=logs/phase2_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Activate venv (run `pip install pandas numpy biopython` if not installed)
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

# Directories
TSV_DIR="../filtered_blast_results"
FASTA_DIR="../cleaned_fastas"
MASTER_COV_DIR="coverage_files"
mkdir -p "$MASTER_COV_DIR"

echo "= Phase 2: Direct capped coverage analysis (0-100%) from filtered TSVs and FASTAs ="

# Embedded Python code for calculation
python_code=$(cat << 'EOF'
import pandas as pd
import numpy as np
import re
from pathlib import Path
from typing import List, Tuple, Dict
from Bio import SeqIO
import sys

def canon_id(s: str) -> str:
    return re.sub(r'\[\]', '', str(s)).split()[0]

def calculate_capped_coverage_pct(intervals: List[Tuple[int, int]], query_length: int) -> float:
    if query_length <= 0 or not intervals:
        return 0.0
    coverage_array = np.zeros(query_length, dtype=bool)
    for start, end in intervals:
        start_idx = max(0, start - 1)
        end_idx = min(query_length, end)
        coverage_array[start_idx:end_idx] = True
    coverage_bases = np.count_nonzero(coverage_array)
    coverage_pct = (coverage_bases / query_length) * 100
    return min(coverage_pct, 100.0)

tsv_path = Path(sys.argv[1])
sample = tsv_path.stem.split('_vs_')[0]
db = tsv_path.stem.split('_vs_')[1].replace('_filtered', '')
fasta_path = Path(sys.argv[2]) / f"{sample}.clean.fasta"
out_csv = Path(sys.argv[3]) / f"{sample}_vs_{db}_capped_stats.csv"

if not fasta_path.exists():
    print(f"SKIP: Missing FASTA for {sample}: {fasta_path}")
    sys.exit(0)

# Load TSV
blast_df = pd.read_csv(tsv_path, sep='\t', names=[
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
])

# Load FASTA lengths
seq_lengths = {canon_id(r.id): len(r.seq) for r in SeqIO.parse(str(fasta_path), 'fasta')}
print(f"Loaded {len(seq_lengths)} sequences for {sample}.")

if blast_df.empty:
    print(f"WARNING: Empty TSV: {tsv_path}")
    results = [{
        'qseqid': qid,
        'num_alignments': 0,
        'pident_mean': 0.0,
        'evalue_min': float('nan'),
        'bitscore_max': 0.0,
        'query_length': qlen,
        'capped_coverage_pct': 0.0
    } for qid, qlen in seq_lengths.items()]
else:
    blast_df['cleaned_qseqid'] = blast_df['qseqid'].apply(canon_id)
    blast_df['query_length'] = blast_df['cleaned_qseqid'].map(seq_lengths)
    blast_df = blast_df.dropna(subset=['query_length'])

    grouped = blast_df.groupby('cleaned_qseqid')
    results = []
    for qid, group in grouped:
        intervals = list(zip(group['qstart'], group['qend']))
        qlen = group['query_length'].max()
        capped_pct = calculate_capped_coverage_pct(intervals, qlen)
        results.append({
            'qseqid': qid,
            'num_alignments': len(group),
            'pident_mean': group['pident'].mean(),
            'evalue_min': group['evalue'].min(),
            'bitscore_max': group['bitscore'].max(),
            'query_length': qlen,
            'capped_coverage_pct': capped_pct
        })

    # Add unaligned
    aligned_ids = {r['qseqid'] for r in results}
    for qid, qlen in seq_lengths.items():
        if qid not in aligned_ids:
            results.append({
                'qseqid': qid,
                'num_alignments': 0,
                'pident_mean': 0.0,
                'evalue_min': float('nan'),
                'bitscore_max': 0.0,
                'query_length': qlen,
                'capped_coverage_pct': 0.0
            })

pd.DataFrame(results).to_csv(out_csv, index=False)
print(f"Saved to {out_csv}")
EOF
)

# Process all TSVs
shopt -s nullglob
for tsv_file in "${TSV_DIR}"/*_filtered.tsv; do
    python3 -c "$python_code" "$tsv_file" "$FASTA_DIR" "$MASTER_COV_DIR"
done

echo "=== Phase 2 complete (outputs in ${MASTER_COV_DIR}) ==="
mkdir -p .checkpoints
: > .checkpoints/phase2_coverage.done
