#!/bin/bash
#SBATCH --job-name=phase2b
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4021
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase2b_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Require seqkit for FASTA subsetting
if ! command -v seqkit >/dev/null 2>&1; then
  echo "ERROR: seqkit not found. Install it via 'conda install -c bioconda seqkit' or ensure it's in PATH." >&2
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$THREADS"
echo "Using THREADS=$THREADS"

# Directories
HIGHCONF_DIR="high_confidence_hits"   # expects ${base}.fna from phase1
OUTPUT_DIR="coverage_filtered"
COVERAGE_ROOT="coverage_files"        # expects subdirs coverage_*
mkdir -p "$OUTPUT_DIR"

# Thresholds
MIN_IDENTITY=0
MIN_COVERAGE=70

# Find coverage stats CSVs inside each coverage_* directory (fixed glob)
shopt -s nullglob
inputs=( "${COVERAGE_ROOT}"/coverage_*/*_coverage_stats.csv )
if (( ${#inputs[@]} == 0 )); then
  echo "ERROR: No *_coverage_stats.csv found under ${COVERAGE_ROOT}/coverage_*" >&2
  exit 1
fi

for stats_csv in "${inputs[@]}"; do
  # stats_csv like: coverage_files/coverage_ES01RB.unmapped_1_classified_bact_90pct/coverage_ES01RB.unmapped_1_classified_bact_90pct_coverage_stats.csv
  file_name="$(basename "$stats_csv" _coverage_stats.csv)"     # coverage_ES01RB.unmapped_1_classified_bact_90pct
  base="${file_name#coverage_}"                                # ES01RB.unmapped_1_classified_bact_90pct

  ids_file="$OUTPUT_DIR/${base}_pass_ids.txt"
  in_fna="$HIGHCONF_DIR/${base}.fna"
  out_fna="$OUTPUT_DIR/${base}_coverage_filtered.fna"
  out_txt="$OUTPUT_DIR/${base}_coverage_filtered_hits.txt"
  done_flag="$OUTPUT_DIR/${base}.done"

  if [[ -f "$done_flag" ]]; then
    echo "Done flag present for $base; skipping."
    continue
  fi

  if [[ ! -s "$in_fna" ]]; then
    echo "WARNING: Missing/empty FASTA: $in_fna; skipping $base"
    : > "$ids_file"
    : > "$out_fna"
    : > "$out_txt"
    : > "$done_flag"
    continue
  fi

  # 1) Build pass ID list from stats
  if [[ -s "$ids_file" ]]; then
    echo "IDs exist for $base, skipping step 1."
  else
    python3 - <<EOF
import pandas as pd, re, os
stats_file = "$stats_csv"
out_final = "$ids_file"
tmp = out_final + ".tmp"

def canon(s: str) -> str:
    return re.sub(r'\\[\\]', '', str(s)).split()[0]

df = pd.read_csv(stats_file)
if df.empty:
    open(tmp, "w").close()
else:
    # Filter by thresholds on per-contig aggregates
    keep = (df['pident_mean'] >= $MIN_IDENTITY) & (df['query_coverage_mean'] >= $MIN_COVERAGE)
    filtered = df.loc[keep].copy()
    filtered['qseqid'] = filtered['qseqid'].map(canon)
    with open(tmp, "w") as f:
        for qid in filtered['qseqid'].dropna().unique():
            f.write(f"{qid}\\n")
os.replace(tmp, out_final)
EOF
  fi

  # If no passing IDs, write empty outputs and continue
  if [[ ! -s "$ids_file" ]]; then
    echo "No sequences passed thresholds for $base; creating empty outputs."
    : > "$out_fna"
    : > "$out_txt"
    : > "$done_flag"
    continue
  fi

  # 2) Subset FASTA with seqkit (normalize headers to canonical IDs before grep)
  if [[ -s "$out_fna" ]]; then
    echo "Filtered FASTA exists for $base, skipping step 2."
  else
    tmp_norm="${in_fna}.norm"
    # Keep only first token of header and strip [] so IDs match the ids_file
    seqkit replace -p '^(\S+).*$' -r '$1' "$in_fna" | sed 's/\[\]//' > "$tmp_norm"
    # Grep using pass IDs (case-insensitive)
    seqkit grep -i -n -f "$ids_file" "$tmp_norm" > "${out_fna}.tmp" || true
    mv "${out_fna}.tmp" "$out_fna"
    rm -f "$tmp_norm"
  fi

  # 3) Summary TSV of passing rows (subset of stats)
  if [[ -s "$out_txt" ]]; then
    echo "Summary exists for $base, skipping step 3."
  else
    python3 - <<EOF2
import pandas as pd, re, os
stats_file = "$stats_csv"
ids_file = "$ids_file"
out_final = "$out_txt"
tmp = out_final + ".tmp"

def canon(s: str) -> str:
    return re.sub(r'\\[\\]', '', str(s)).split()[0]

# Read pass IDs
with open(ids_file) as f:
    pass_ids = set(line.strip() for line in f if line.strip())

df = pd.read_csv(stats_file)
if df.empty or not pass_ids:
    df_out = df.head(0)
else:
    df['qseqid'] = df['qseqid'].map(canon)
    df_out = df[df['qseqid'].isin(pass_ids)]

df_out.to_csv(tmp, sep="\\t", index=False)
os.replace(tmp, out_final)
EOF2
  fi

  passed_count=$(wc -l < "$ids_file" | awk '{print $1}')
  echo "Filtered $passed_count sequences; outputs written for $base"
  : > "$done_flag"
done


