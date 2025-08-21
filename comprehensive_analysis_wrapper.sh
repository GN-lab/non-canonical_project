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

# Check for coverage data (essential)
if [[ ! -f "cov_analysis/coverage_depth_stats.tsv" ]]; then
    echo "ERROR: Coverage file not found: cov_analysis/coverage_depth_stats.tsv" >&2
    echo "Please run Phase 4 first to generate coverage statistics." >&2
    exit 1
fi

# Check for functional annotations (optional but nice to have)
if ! ls functional_analysis/*.tsv >/dev/null 2>&1; then
    echo "WARNING: No functional annotation files found in functional_analysis/"
    echo "Phase 5 may have failed or not been run yet."
    echo "Continuing with coverage analysis only..."
fi

echo "Starting Phase 6: Comprehensive Analysis"

# Run Python analysis (primary)
echo "Running Python comprehensive analysis..."
python comprehensive_analysis.py

# Try to run R analysis if available (but don't fail if it doesn't work)
if [[ -f "comprehensive_analysis.r" ]]; then
    echo "Running R analysis (if available)..."
    # Use timeout to prevent hanging and continue even if R fails
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
    
    # Print summary
    echo ""
    echo "=== PHASE 6 COMPLETE ==="
    echo "Python outputs: comprehensive_plots_python/"
    if [[ -d "comprehensive_plots" ]]; then
        echo "R outputs: comprehensive_plots/"
    fi
    echo "Checkpoint: .checkpoints/phase6_comprehensive.done"
else
    echo "ERROR: No analysis outputs found!" >&2
    echo "Please check the log files for errors." >&2
    exit 1
fi
