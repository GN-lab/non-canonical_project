#!/bin/bash
#SBATCH --job-name=phase3
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8042
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase3_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Activate the environment that has your Python deps (choose one)
# conda activate base
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

# Directories
COV_FILTERED_DIR="coverage_filtered"  # contains *_coverage_filtered_hits.txt
KRAKEN_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/kraken_results"  # *_kraken_output_with_name.txt
OUTPUT_DIR="intersection_results"

# Checks
[[ -d "$COV_FILTERED_DIR" ]] || { echo "ERROR: Coverage filtered dir missing: $COV_FILTERED_DIR" >&2; exit 1; }
[[ -d "$KRAKEN_DIR" ]] || { echo "ERROR: Kraken dir missing: $KRAKEN_DIR" >&2; exit 1; }
if ! ls "$COV_FILTERED_DIR"/*_coverage_filtered_hits.txt >/dev/null 2>&1; then
  echo "ERROR: No filtered hits TXT in $COV_FILTERED_DIR; Phase 2b must have failed." >&2
  exit 1
fi

echo "Intersecting using: BLAST hits dir=$COV_FILTERED_DIR, Kraken dir=$KRAKEN_DIR, Output dir=$OUTPUT_DIR"

# Optional: export COMBINE=1 to create a single TSV
COMBINE="${COMBINE:-0}"
if [[ "$COMBINE" -eq 1 ]]; then
  COMBINE_FLAG="--combine"
else
  COMBINE_FLAG=""
fi

python3 intersect_kraken_blast.py \
  --blast_hits_dir "$COV_FILTERED_DIR" \
  --kraken_dir "$KRAKEN_DIR" \
  --output_dir "$OUTPUT_DIR" \
  $COMBINE_FLAG

if ls "$OUTPUT_DIR"/*_kraken_blast_intersect.tsv >/dev/null 2>&1 || { [[ "$COMBINE" -eq 1 ]] && [[ -s "$OUTPUT_DIR/all_intersections.tsv" ]]; }; then
  mkdir -p .checkpoints
  : > .checkpoints/phase3_intersect.done
  echo "Phase 3 checkpoint written."
else
  echo "WARNING: No Phase 3 outputs detected in $OUTPUT_DIR; not writing checkpoint."
fi
