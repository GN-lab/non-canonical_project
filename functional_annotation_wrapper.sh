#!/bin/bash
#SBATCH --job-name=phase5
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8042
#SBATCH --time=24:00:00
#SBATCH --output=logs/phase5_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

source "/data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh"

INPUT_DIR="length_filtered"
# Require Phase 4 outputs
if ! ls "$INPUT_DIR"/*.fna >/dev/null 2>&1; then
  echo "ERROR: No FASTA inputs in $INPUT_DIR; Phase 4 must have failed." >&2
  exit 1
fi

# Run Phase 5
bash functional_annotation.sh

# Only write checkpoint if expected Phase 5 outputs exist
# Adjust these patterns to match your script's actual outputs
if ls functional_analysis/*_nr.tsv >/dev/null 2>&1 ||
   ls functional_analysis/*_pfam.tbl >/dev/null 2>&1 ||
   ls functional_analysis/*_cog.tsv >/dev/null 2>&1; then
  mkdir -p .checkpoints
  : > .checkpoints/phase5_functional.done
  echo "Phase 5 checkpoint written."
else
  echo "WARNING: No Phase 5 outputs detected; not writing checkpoint."
fi

