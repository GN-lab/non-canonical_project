#!/usr/bin/env bash
set -euo pipefail

BASE='../../nextflow_rnaseq/output_hartwig/star_salmon/unmapped/'
OUTDIR='merged_fasta'
mkdir -p "$OUTDIR"

shopt -s nullglob
r1_files=("$BASE"/*.unmapped_1.fastq)
if (( ${#r1_files[@]} == 0 )); then
  echo "No R1 files found under $BASE (*.unmapped_1.fastq)" >&2
  exit 1
fi

for r1 in "${r1_files[@]}"; do
  sample=$(basename "$r1" .unmapped_1.fastq)
  r2="$BASE/$sample.unmapped_2.fastq"

  if [[ ! -f "$r2" ]]; then
    echo "Skipping $sample: missing $r2" >&2
    continue
  fi

  r2rc="$OUTDIR/${sample}_R2_rc.fastq"
  merged="$OUTDIR/${sample}_merged.fastq"  # Change to .fastq.gz if you want compressed output

  echo "[$sample] Reverse complementing R2 -> $r2rc ..."
  # seqkit detects gz input and writes gz when -o ends with .gz
  seqkit seq -r -p -t DNA "$r2" -o "$r2rc"

  echo "[$sample] Merging R1 and R2_rc -> $merged ..."
  # Use cat instead of zcat since files are uncompressed; optionally gzip the output
  cat -- "$r1" "$r2rc" > "$merged"  # For uncompressed output
  # Or: cat -- "$r1" "$r2rc" | gzip -c > "${merged}.gz"  # For compressed output (rename accordingly)

  echo "[$sample] Read counts verification:"
  # Replace zcat with wc directly for uncompressed files
  r1_n=$(( $(wc -l < "$r1") / 4 ))
  r2_n=$(( $(wc -l < "$r2") / 4 ))
  r2rc_n=$(( $(wc -l < "$r2rc") / 4 ))
  merged_n=$(( $(wc -l < "$merged") / 4 ))  # Adjust if you gzipped the merged file

  echo "  Original R1: $r1_n"
  echo "  Original R2: $r2_n"
  echo "  R2 reverse complement: $r2rc_n"
  echo "  Merged total: $merged_n"

  if (( merged_n != r1_n + r2_n )); then
    echo "  WARNING: merged count ($merged_n) != R1+R2 ($((r1_n + r2_n)))" >&2
  fi
done
