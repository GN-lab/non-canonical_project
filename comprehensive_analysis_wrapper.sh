#!/bin/bash
#SBATCH --job-name=phase6
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4  # Increased for better performance
#SBATCH --mem-per-cpu=8042  # Increased for large datasets
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase6_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Load modules and activate environment
module load R
export R_LIBS_USER="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/R_libs"
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

# Check for required files
echo "Checking for required files..."

# Check for coverage data (essential; prefer filtered if available)
COVERAGE_FILE="cov_analysis/coverage_depth_stats.tsv"
FILTERED_DIR="coverage_filtered"
if [[ -d "$FILTERED_DIR" && $(ls "$FILTERED_DIR"/*_coverage_filtered_hits.txt 2>/dev/null | wc -l) -gt 0 ]]; then
    echo "Using filtered coverage data from $FILTERED_DIR"
    # Compile filtered stats if needed (optional; assume already compiled or add compilation here)
    COVERAGE_FILE="cov_analysis/filtered_coverage_stats.tsv"  # Adjust if you have a compiled file
else
    if [[ ! -f "$COVERAGE_FILE" ]]; then
        echo "ERROR: Coverage file not found: $COVERAGE_FILE" >&2
        echo "Please run Phase 4 or 2b first to generate coverage statistics." >&2
        exit 1
    fi
    echo "Using raw coverage data: $COVERAGE_FILE"
fi

# Check for functional annotations (optional)
if ! ls functional_analysis/*.tsv >/dev/null 2>&1; then
    echo "WARNING: No functional annotation files found in functional_analysis/"
    echo "Continuing with coverage analysis only..."
fi

echo "Starting Phase 6: Comprehensive Analysis"

# Run Python analysis (primary; updated with E-value integration below)
echo "Running Python comprehensive analysis..."
python comprehensive_analysis.py

# Try to run R analysis if available (but don't fail if it doesn't work)
if [[ -f "comprehensive_analysis.r" ]]; then
    echo "Running R analysis (if available)..."
    timeout 1h Rscript comprehensive_analysis.r || {
        echo "WARNING: R analysis failed or timed out. Continuing with Python results..."
    }
else
    echo "No R script found. Using Python analysis only."
fi

# Check for outputs and write checkpoint
echo "Checking for analysis outputs..."

# Check for Python outputs
python_outputs=(
    "comprehensive_plots_python/scatter_alignments_main.png"
    "comprehensive_plots_python/scatter_alignments_by_group.png"
    "comprehensive_plots_python/alignment_statistics_detailed.csv"
)

output_found=false
for output in "${python_outputs[@]}"; do
    if [[ -f "$output" ]]; then
        echo "Found output: $output"
        output_found=true
    fi
done

# Also check for R outputs if they exist
if [[ -d "comprehensive_plots" ]]; then
    echo "Found R outputs directory: comprehensive_plots/"
    output_found=true
fi

if [[ "$output_found" == true ]]; then
    mkdir -p .checkpoints
    : > .checkpoints/phase6_comprehensive.done
    echo "Phase 6 checkpoint written."
else
    echo "ERROR: No analysis outputs found!" >&2
    echo "Please check the log files for errors." >&2
    exit 1
fi
