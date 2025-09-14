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

###############################################################################
# Requirements
###############################################################################
if ! command -v seqkit >/dev/null 2>&1; then
  echo "ERROR: seqkit not found.  Install with: conda install -c bioconda seqkit" >&2
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$THREADS"
echo "Using THREADS=$THREADS"

###############################################################################
# Paths and thresholds
###############################################################################
FASTA_DIR="../cleaned_fastas"     # original FASTAs ( one SAMPLE.clean.fasta each )
COVERAGE_ROOT="coverage_files"    # holds coverage_*/ *_capped_stats.csv
OUTPUT_DIR="coverage_filtered"
mkdir -p "$OUTPUT_DIR"

# Thresholds
MIN_IDENTITY=80          # set to 80-90 if you want close matches
MIN_COVERAGE=70         # % query coverage (capped_coverage_pct)
MIN_EVALUE=1e-5         # maximum allowed e-value

###############################################################################
# Locate input CSVs
###############################################################################
shopt -s nullglob
inputs=( "${COVERAGE_ROOT}"/*_capped_stats.csv )
if (( ${#inputs[@]} == 0 )); then
  echo "ERROR: No *_capped_stats.csv found under ${COVERAGE_ROOT}/coverage_*" >&2
  exit 1
fi

###############################################################################
# Main loop
###############################################################################
for stats_csv in "${inputs[@]}"; do
  # e.g. coverage_files/coverage_HT01RB_bact/coverage_HT01RB_bact_capped_stats.csv
  file_name=$(basename "$stats_csv" _capped_stats.csv)   # coverage_HT01RB_bact
  base=${file_name#coverage_}                            # HT01RB_bact
  sample=${base%%_*}                                     # HT01RB  (assumes base=SAMPLE_db)

  ids_file="$OUTPUT_DIR/${base}_pass_ids.txt"
  in_fna="$FASTA_DIR/${sample}.clean.fasta"
  out_fna="$OUTPUT_DIR/${base}_coverage_filtered.fna"
  out_txt="$OUTPUT_DIR/${base}_coverage_filtered_hits.txt"
  done_flag="$OUTPUT_DIR/${base}.done"

  [[ -f "$done_flag" ]] && { echo "Already done: $base"; continue; }

  if [[ ! -s "$in_fna" ]]; then
    echo "WARNING: FASTA missing for $sample: $in_fna"
    : > "$ids_file" "$out_fna" "$out_txt" "$done_flag"
    continue
  fi

  ###########################################################################
  # 1) create list of qseqids that pass thresholds
  ###########################################################################
  if [[ ! -s "$ids_file" ]]; then
    python3 - <<EOF
import pandas as pd, re, os
stats_file = "$stats_csv"
tmp_out = "$ids_file.tmp"

canon = lambda s: re.sub(r'\[\]', '', str(s)).split()[0]

df = pd.read_csv(stats_file)
if not df.empty:
    keep = (
        (df['pident_mean'] >= $MIN_IDENTITY) &
        (df['capped_coverage_pct'] >= $MIN_COVERAGE) &
        (df['evalue_min'] <= $MIN_EVALUE)
    )
    df = df.loc[keep, 'qseqid'].map(canon).dropna().unique()
else:
    df = []

with open(tmp_out, 'w') as fh:
    for qid in df:
        fh.write(f"{qid}\n")

os.replace(tmp_out, "$ids_file")
EOF
  fi

  # nothing passed?  make empty outputs and continue
  if [[ ! -s "$ids_file" ]]; then
    echo "No sequences passed thresholds for $base"
    : > "$out_fna" "$out_txt" "$done_flag"
    continue
  fi

  ###########################################################################
  # 2) subset FASTA with seqkit
  ###########################################################################
  if [[ ! -s "$out_fna" ]]; then
    tmp_norm="${in_fna}.norm"
    # keep only first token in header -> matches qseqid
    seqkit replace -p '^(\S+).*$' -r '$1' "$in_fna" | sed 's/\[\]//' > "$tmp_norm"
    seqkit grep -i -n -f "$ids_file" "$tmp_norm" > "${out_fna}.tmp" || true
    mv "${out_fna}.tmp" "$out_fna"
    rm -f "$tmp_norm"
  fi

  ###########################################################################
  # 3) write summary TSV of passing rows
  ###########################################################################
  if [[ ! -s "$out_txt" ]]; then
    python3 - <<EOF
import pandas as pd, re, os, sys
stats_file = "$stats_csv"
ids      = {line.strip() for line in open("$ids_file") if line.strip()}
tmp_out  = "$out_txt.tmp"

canon = lambda s: re.sub(r'\[\]', '', str(s)).split()[0]
df = pd.read_csv(stats_file)
if not df.empty and ids:
    df['qseqid'] = df['qseqid'].map(canon)
    df = df[df['qseqid'].isin(ids)]
else:
    df = df.head(0)   # empty with header

df.to_csv(tmp_out, sep="\t", index=False)
os.replace(tmp_out, "$out_txt")
EOF
  fi

  echo "Finished $base ($(wc -l < "$ids_file") sequences kept)"
  : > "$done_flag"
done
