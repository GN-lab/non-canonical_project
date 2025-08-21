#!/bin/bash
#SBATCH --job-name=phase4
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4021
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase4_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

module load R
export R_LIBS_USER="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/R_libs"
# Activate the virtual environment (ensure biopython, pandas/numpy are available)
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

# Require Phase 2b outputs (check if coverage_filtered dir has files)
COVERAGE_DIR="coverage_filtered"
if ! ls "$COVERAGE_DIR"/*_coverage_filtered.fna >/dev/null 2>&1; then
  echo "ERROR: No *_coverage_filtered.fna files found in $COVERAGE_DIR." >&2
  exit 1
fi

# Run Phase 4 on Phase 2b files
bash length_filter.sh

# Only checkpoint if Phase 4 outputs exist
if ls cov_analysis/*_hits.tsv >/dev/null 2>&1; then
  mkdir -p .checkpoints
  : > .checkpoints/phase4_length_filter.done
  echo "Phase 4 checkpoint written."
else
  echo "WARNING: No Phase 4 outputs detected; not writing checkpoint."
fi
