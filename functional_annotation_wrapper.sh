#!/bin/bash
#SBATCH --job-name=phase5
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8042
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/phase5_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Optional: activate conda profile if needed for environment setup
source "/data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh"

# Configuration
INPUT_DIR="length_filtered"
OUTPUT_DIR="functional_analysis"
WORKER="functional_annotation.sh"
LOCK_FILE=".locks/phase5.lock"

# Ensure directories
mkdir -p "$OUTPUT_DIR" .locks logs .checkpoints

# Verify Phase 4 outputs exist
if ! ls "$INPUT_DIR"/*.fna >/dev/null 2>&1; then
  echo "ERROR: No FASTA inputs in $INPUT_DIR; Phase 4 must have failed." >&2
  exit 1
fi

# Acquire a non-blocking lock to prevent overlapping runs
exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  echo "Another phase5 wrapper instance is running; exiting."
  exit 0
fi

echo "Scanning for pending samples..."

# Build the list of samples that are incomplete (missing or empty outputs)
shopt -s nullglob
pending=()
for fasta in "$INPUT_DIR"/*.fna; do
  [[ -s "$fasta" ]] || continue
  sample=$(basename "$fasta" .fna)
  prefix="${OUTPUT_DIR}/${sample}"

  # Expected outputs indicating a complete sample
  expected=(
    "${prefix}_nr.tsv"
    "${prefix}_translated.faa"
    "${prefix}_pfam.tbl"
    "${prefix}_pfam_domains.tbl"
    "${prefix}_cog.tsv"
  )

  complete=true
  for f in "${expected[@]}"; do
    if [[ ! -s "$f" ]]; then
      complete=false
      break
    fi
  done

  if ! $complete; then
    pending+=("$fasta")
  fi
done

if (( ${#pending[@]} == 0 )); then
  echo "All samples complete; nothing to do."
else
  echo "Pending samples: ${#pending[@]}"

  # Process each pending sample by isolating it in a temporary INPUT_DIR
  for fasta in "${pending[@]}"; do
    sample=$(basename "$fasta" .fna)
    echo "Processing sample: $sample"

    tmpdir=$(mktemp -d)
    # Ensure cleanup on error within this loop iteration
    cleanup() { rm -rf "$tmpdir"; }
    trap cleanup RETURN

    cp -f "$fasta" "$tmpdir"/

    # Run the worker with the temp directory so it only sees this one sample
    INPUT_DIR="$tmpdir" bash "$WORKER"

    # Cleanup temp
    rm -rf "$tmpdir"
    trap - RETURN
  done
fi

# Write checkpoint only if there are non-empty outputs
if find "$OUTPUT_DIR" -type f \( -name "*_nr.tsv" -o -name "*_pfam.tbl" -o -name "*_cog.tsv" \) ! -size 0 | grep -q .; then
  : > .checkpoints/phase5_functional.done
  echo "Phase 5 checkpoint written."
else
  echo "WARNING: No complete Phase 5 outputs detected; not writing checkpoint."
fi
