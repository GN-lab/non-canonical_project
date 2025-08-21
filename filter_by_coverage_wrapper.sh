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

source "/data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh"
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

shopt -s nullglob

# Corrected glob (no stray underscore)
stats=( coverage_files/coverage_*/*_coverage_stats.csv )

# Existence check in the new location
if ! ls coverage_files/coverage_*/*_coverage_stats.csv >/dev/null 2>&1; then
  echo "ERROR: No *_coverage_stats.csv (Phase 2 outputs) found under coverage_files" >&2
  exit 1
fi

# Run the filter
bash filter_by_coverage.sh

# Checkpoint: look for outputs where Phase 2b writes them (coverage_filtered)
if ls coverage_filtered/*_coverage_filtered.csv >/dev/null 2>&1 || ls coverage_filtered/*_kept.txt >/dev/null 2>&1; then
  mkdir -p .checkpoints
  : > .checkpoints/phase2b_filter_by_cov.done
  echo "Phase 2b checkpoint written."
else
  echo "WARNING: No Phase 2b outputs detected; not writing checkpoint."
fi
