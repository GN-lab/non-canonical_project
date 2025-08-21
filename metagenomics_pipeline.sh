#!/bin/bash
#SBATCH --job-name=metagenomics_driver
#SBATCH --partition=master-worker
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4021
#SBATCH --time=48:00:00
#SBATCH --output=logs/master_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "**ERROR** rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$THREADS"

# Shared environment
module load R
export R_LIBS_USER="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/R_libs"
source /data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh

mkdir -p logs .checkpoints
export BLAST_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results"

echo "=== Submitting Metagenomic Pipeline (dependency-driven) ==="

[[ -d "$BLAST_DIR" ]] || { echo "ERROR: BLAST_DIR missing: $BLAST_DIR" >&2; exit 1; }

# Phase completeness tests: prefer real outputs, fall back to checkpoints
phase1_done() { ls high_confidence_hits/*_90pct.tsv >/dev/null 2>&1 || [[ -f .checkpoints/phase1_extract.done ]]; }
phase2_done() { ls coverage_*/*_coverage_stats.csv >/dev/null 2>&1 || [[ -f .checkpoints/phase2_coverage.done ]]; }
phase2b_done(){ ls coverage_filtered/*_coverage_filtered.fna >/dev/null 2>&1 || [[ -f .checkpoints/phase2b_filter_by_cov.done ]]; }
phase3_done() { ls joint_classified* >/dev/null 2>&1 || [[ -f .checkpoints/phase3_intersect.done ]]; }
phase4_done() { ls length_filtered/*.fna >/dev/null 2>&1 || [[ -f .checkpoints/phase4_length_filter.done ]]; }
phase5_done() { ls functional_analysis/*_nr.tsv >/dev/null 2>&1 || [[ -f .checkpoints/phase5_functional.done ]]; }
phase6_done() { ls results/* >/dev/null 2>&1 || ls reports/* >/dev/null 2>&1 || [[ -f .checkpoints/phase6_comprehensive.done ]]; }

# Rolling dependency: the last job submitted in THIS run
jobPrev=""

# ---- Phase 1 ----
job1=""
if phase1_done; then
  echo "Phase 1 already submitted/completed."
else
  echo "Phase 1: Extracting 90% identity hits..."
  # Sanity: confirm inputs exist
  if ! ls "$BLAST_DIR"/classified/*_classified_vs_*.tsv >/dev/null 2>&1 && \
     ! ls "$BLAST_DIR"/unclassified/*_unclassified_vs_*.tsv >/dev/null 2>&1; then
    echo "ERROR: No BLAST TSVs in $BLAST_DIR/classified or unclassified" >&2
    exit 1
  fi
  job1=$(sbatch --parsable extract_90percent_hits.sh)
  echo "Submitted Phase 1 as $job1"
  jobPrev="$job1"
fi

# ---- Phase 2 ----
job2=""
if phase2_done; then
  echo "Phase 2 already submitted/completed."
else
  echo "Phase 2: Coverage analysis..."
  if [[ -n "${jobPrev:-}" ]]; then
    job2=$(sbatch --parsable --dependency=afterok:${jobPrev} coverage_phase2_wrapper.sh)
    echo "Submitted Phase 2 as $job2 (afterok:$jobPrev)"
  else
    job2=$(sbatch --parsable coverage_phase2_wrapper.sh)
    echo "Submitted Phase 2 as $job2"
  fi
  jobPrev="$job2"
fi

# ---- Phase 2b ----
job2b=""
if phase2b_done; then
  echo "Phase 2b already submitted/completed."
else
  echo "Phase 2b: Coverage filtering..."
  if [[ -n "${jobPrev:-}" ]]; then
    job2b=$(sbatch --parsable --dependency=afterok:${jobPrev} filter_by_coverage_wrapper.sh)
    echo "Submitted Phase 2b as $job2b (afterok:$jobPrev)"
  else
    job2b=$(sbatch --parsable filter_by_coverage_wrapper.sh)
    echo "Submitted Phase 2b as $job2b"
  fi
  jobPrev="$job2b"
fi

# ---- Phase 3 ----
job3=""
if phase3_done; then
  echo "Phase 3 already submitted/completed."
else
  echo "Phase 3: Kraken2 ∩ BLAST intersect..."
  if [[ -n "${jobPrev:-}" ]]; then
    job3=$(sbatch --parsable --dependency=afterok:${jobPrev} intersect_phase3_wrapper.sh)
    echo "Submitted Phase 3 as $job3 (afterok:$jobPrev)"
  else
    job3=$(sbatch --parsable intersect_phase3_wrapper.sh)
    echo "Submitted Phase 3 as $job3"
  fi
  jobPrev="$job3"
fi

# ---- Phase 4 ----
job4=""
if phase4_done; then
  echo "Phase 4 already submitted/completed."
else
  echo "Phase 4: Length filtering..."
  if [[ -n "${jobPrev:-}" ]]; then
    job4=$(sbatch --parsable --dependency=afterok:${jobPrev} length_filter_wrapper.sh)
    echo "Submitted Phase 4 as $job4 (afterok:$jobPrev)"
  else
    job4=$(sbatch --parsable length_filter_wrapper.sh)
    echo "Submitted Phase 4 as $job4"
  fi
  jobPrev="$job4"
fi

# ---- Phase 5 ----
job5=""
if phase5_done; then
  echo "Phase 5 already submitted/completed."
else
  echo "Phase 5: Functional annotation..."
  if [[ -n "${jobPrev:-}" ]]; then
    job5=$(sbatch --parsable --dependency=afterok:${jobPrev} functional_annotation_wrapper.sh)
    echo "Submitted Phase 5 as $job5 (afterok:$jobPrev)"
  else
    job5=$(sbatch --parsable functional_annotation_wrapper.sh)
    echo "Submitted Phase 5 as $job5"
  fi
  jobPrev="$job5"
fi

# ---- Phase 6 ----
job6=""
if phase6_done; then
  echo "Phase 6 already submitted/completed."
else
  echo "Phase 6: Comprehensive analysis in R..."
  if [[ -n "${jobPrev:-}" ]]; then
    job6=$(sbatch --parsable --dependency=afterok:${jobPrev} comprehensive_analysis_wrapper.sh)
    echo "Submitted Phase 6 as $job6 (afterok:$jobPrev)"
  else
    job6=$(sbatch --parsable comprehensive_analysis_wrapper.sh)
    echo "Submitted Phase 6 as $job6"
  fi
  jobPrev="$job6"
fi

echo "=== Submission complete. Jobs chained with afterok dependencies (where applicable). ==="
