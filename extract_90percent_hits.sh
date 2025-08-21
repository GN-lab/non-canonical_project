#!/bin/bash
#SBATCH --job-name=phase1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8042
#SBATCH --time=24:00:00
#SBATCH --output=logs/phase1_%j.log

set -eo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$THREADS"
echo "Using THREADS=$THREADS"

# Adjust these to your layout
BLAST_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results"
CLASSIFIED_TSV_DIR="${BLAST_DIR}/classified"
UNCLASSIFIED_TSV_DIR="${BLAST_DIR}/unclassified"

# NEW: Directories for classified/unclassified FASTA files (query inputs to BLAST)
CLASSIFIED_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/kraken_results/classified"
UNCLASSIFIED_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/kraken_results/unclassified"

OUTPUT_DIR="high_confidence_hits"
mkdir -p "$OUTPUT_DIR"

atomic_mv() {
  local tmp="$1" final="$2"
  mv -f "$tmp" "$final"
}

filter_tsv() {
  local sample="$1"   # e.g., ES01RB.unmapped_1
  local type="$2"     # classified | unclassified
  local db="$3"       # bact | fungi | parasite | virus | etc.
  local tsv="$4"

  echo "Processing sample='${sample}' type='${type}' db='${db}'..."

  if [[ ! -s "$tsv" ]]; then
    echo "  WARNING: Missing/empty BLAST TSV: $tsv; skipping."
    return
  fi

  # Locate corresponding FASTA file based on type (assume named like ${sample}.fasta)
  local fasta_dir
  if [[ "$type" == "classified" ]]; then
    fasta_dir="$CLASSIFIED_DIR"
  elif [[ "$type" == "unclassified" ]]; then
    fasta_dir="$UNCLASSIFIED_DIR"
  else
    echo "  ERROR: Unknown type '${type}'; skipping."
    return
  fi
  local fasta="${fasta_dir}/${sample}.fasta"
  if [[ ! -s "$fasta" ]]; then
    echo "  ERROR: Missing/empty FASTA for sample '${sample}' in ${fasta_dir}: $fasta; skipping."
    return
  fi

  # Compute query_lengths.txt from FASTA (qseqid\tlength)
  local query_len_file="${OUTPUT_DIR}/${sample}_${type}_${db}_query_lengths.txt"
  if [[ ! -s "$query_len_file" ]]; then
    awk '/^>/ {if (id) print id "\t" len; id = substr($1, 2); len=0; next} {len += length($0)} END {if (id) print id "\t" len}' "$fasta" > "$query_len_file"
    echo "  Generated query_lengths.txt from FASTA."
  fi

  local hits_tsv="${OUTPUT_DIR}/${sample}_${type}_${db}_90pct.tsv"
  local qids="${OUTPUT_DIR}/${sample}_${type}_${db}_qids.txt"
  # CHANGED: write sequences to _90pct.fna so downstream expects real nucleotide sequences there
  local out_fna="${OUTPUT_DIR}/${sample}_${type}_${db}_90pct.fna"

  if [[ -s "$hits_tsv" ]]; then
    echo "  Filtered hits exist, skipping filtering."
  else
    local tmp_hits="${hits_tsv}.tmp"
    # Filters: min query length >=20; pident>=80; aln length>=50; mismatches<=5; bitscore>=50
    awk -F'\t' -v len_file="$query_len_file" '
      BEGIN {
        while ((getline < len_file) > 0) {
          query_len[$1] = $2
        }
        close(len_file)
      }
      $1 in query_len && query_len[$1] >= 20 &&
      $3 >= 80 && $4 >= 50 && $5 <= 5 && $12 >= 50 {
        print $0 "\t" query_len[$1]  # Append query_length as new column
      }
    ' "$tsv" > "$tmp_hits"
    atomic_mv "$tmp_hits" "$hits_tsv"
  fi

  if [[ ! -s "$hits_tsv" ]]; then
    : > "$qids"
    echo "  No high-confidence hits in $(basename "$tsv")."
    return
  fi

  if [[ -s "$qids" ]]; then
    echo "  qids already exist, skipping qid extraction."
  else
    local tmp_qids="${qids}.tmp"
    cut -f1 "$hits_tsv" | sort -u > "$tmp_qids"
    atomic_mv "$tmp_qids" "$qids"
  fi

  # Extract nucleotide sequences to .fna using qids and original fasta
  if [[ -s "$out_fna" ]]; then
    echo "  .fna exists, skipping extraction."
  else
    local tmp_fna="${out_fna}.tmp"
    # Inline Python to extract sequences matching qids (canonicalized)
    python3 - "$fasta" "$qids" "$tmp_fna" <<'PY'
import sys
from Bio import SeqIO

fasta = sys.argv[1]
qids_file = sys.argv[2]
out_fna = sys.argv[3]

# Load qids (canonical: strip [] and after space)
qids = set()
with open(qids_file, 'r') as f:
    for line in f:
        id = line.strip().split()[0].replace('[]', '')
        qids.add(id)

# Parse FASTA and write matching records
with open(out_fna, 'w') as out:
    for record in SeqIO.parse(fasta, 'fasta'):
        rec_id = record.id.split()[0].replace('[]', '')
        if rec_id in qids:
            SeqIO.write(record, out, 'fasta')
PY
    atomic_mv "$tmp_fna" "$out_fna"
    echo "  Extracted nucleotide sequences to $(basename "$out_fna")"
  fi

  echo "  → Filtered: $(basename "$hits_tsv"); qids: $(basename "$qids"); fna: $(basename "$out_fna")"
}

process_file() {
  local tsv="$1"
  local fname; fname="$(basename "$tsv")"

  # Split strictly on the exact tokens
  local left right type
  if [[ "$fname" == *_classified_vs_* ]]; then
    type="classified"
    left="${fname%%_classified_vs_*}"   # before token
    right="${fname#*_classified_vs_}"   # after token
  elif [[ "$fname" == *_unclassified_vs_* ]]; then
    type="unclassified"
    left="${fname%%_unclassified_vs_*}"
    right="${fname#*_unclassified_vs_}"
  else
    echo "WARNING: Unexpected filename: $fname"
    return
  fi

  local sample="$left"
  local db="${right%.tsv}"

  if [[ -z "$sample" || -z "$type" || -z "$db" ]]; then
    echo "WARNING: Failed to parse: $fname -> sample='$sample' type='$type' db='$db'"
    return
  fi

  echo "Parsed -> sample='${sample}' type='${type}' db='${db}' from '${fname}'"
  filter_tsv "$sample" "$type" "$db" "$tsv"
}

main() {
  shopt -s nullglob
  local processed_any=0

  for tsv in "$CLASSIFIED_TSV_DIR"/*_classified_vs_*.tsv "$UNCLASSIFIED_TSV_DIR"/*_unclassified_vs_*.tsv; do
    [[ -f "$tsv" ]] || continue
    process_file "$tsv"
    processed_any=1
  done

  if [[ "$processed_any" -eq 0 ]]; then
    echo "ERROR: No *_classified_vs_*.tsv or *_unclassified_vs_*.tsv found under:"
    echo "  - $CLASSIFIED_TSV_DIR"
    echo "  - $UNCLASSIFIED_TSV_DIR"
    exit 1
  fi

  # Reaching here means the script ran without hitting the "no inputs" error
  mkdir -p .checkpoints
  : > .checkpoints/phase1_extract.done
}

main "$@"
