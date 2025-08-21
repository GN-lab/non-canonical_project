#!/bin/bash
#SBATCH --job-name=phase2
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8042
#SBATCH --time=24:00:00
#SBATCH --output=logs/phase2_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Environment
module load R
export R_LIBS_USER="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/R_libs"

# Activate the virtual environment with matplotlib installed
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

# Absolute source directories
CLASSIFIED_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/kraken_results/classified"
UNCLASSIFIED_DIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/kraken_results/unclassified"

# Phase 1 outputs
HCH_DIR="high_confidence_hits"

# Output directory for assembled .fna files
ASSEMBLED_DIR="assembled_fnas"
mkdir -p "$ASSEMBLED_DIR"

# Master coverage results directory (all per-sample subdirs will live here)
MASTER_COV_DIR="coverage_files"
mkdir -p "$MASTER_COV_DIR"

# Check that seqkit is available
if ! command -v seqkit >/dev/null 2>&1; then
  echo "ERROR: seqkit not found on PATH. Ensure conda environment is activated correctly." >&2
  exit 1
fi

shopt -s nullglob
tsvs=( ${HCH_DIR}/*_90pct.tsv )
if (( ${#tsvs[@]} == 0 )); then
  echo "ERROR: No Phase 1 outputs found: ${HCH_DIR}/*_90pct.tsv" >&2
  exit 1
fi

echo "= Phase 2: Building per-TSV FASTAs from qids + source FASTAs (robust parse, seqkit) ="

trim_suffix() {
  local s="$1" suf="$2"
  if [[ "$s" == *"$suf" ]]; then
    echo "${s%"$suf"}"
  else
    echo "$s"
  fi
}

parse_base() {
  # Input: basename like "<sample>_<type>_<db>_90pct"
  # Output globals: PAR_SAMPLE, PAR_TYPE, PAR_DB
  local base="$1"
  PAR_SAMPLE=""; PAR_TYPE=""; PAR_DB=""

  local no90
  no90="$(trim_suffix "$base" "_90pct")"

  local dbs=("_bact" "_fungi" "_parasite" "_virus")
  local found_db=""
  for d in "${dbs[@]}"; do
    if [[ "$no90" == *"$d" ]]; then
      found_db="$d"
      break
    fi
  done
  if [[ -z "$found_db" ]]; then
    return 1
  fi
  local no_db="${no90%"$found_db"}"
  PAR_DB="${found_db#_}"

  local type=""
  if [[ "$no_db" == *_classified ]]; then
    type="classified"
    PAR_SAMPLE="${no_db%_classified}"
  elif [[ "$no_db" == *_unclassified ]]; then
    type="unclassified"
    PAR_SAMPLE="${no_db%_unclassified}"
  else
    return 1
  fi
  PAR_TYPE="$type"

  [[ -n "$PAR_SAMPLE" && -n "$PAR_TYPE" && -n "$PAR_DB" ]] || return 1
  return 0
}

find_source_fasta() {
  local type="$1" sample="$2"
  local base_dir
  case "$type" in
    classified)   base_dir="$CLASSIFIED_DIR" ;;
    unclassified) base_dir="$UNCLASSIFIED_DIR" ;;
    *) return 1 ;;
  esac

  local cand
  for ext in fasta fa fna; do
    cand="${base_dir}/${sample}.${ext}"
    if [[ -s "$cand" ]]; then
      echo "$cand"
      return 0
    fi
  done
  return 1
}

build_fasta_for_base() {
  local base="$1"  # sample_type_db_90pct
  local qids="${HCH_DIR}/${base%_90pct}_qids.txt"
  local out_fa="${HCH_DIR}/${base}.fna"

  if ! parse_base "$base"; then
    echo "  SKIP: Parse failed for base='$base'"
    return 1
  fi

  local sample="$PAR_SAMPLE"
  local type="$PAR_TYPE"
  local db="$PAR_DB"

  echo "  Build: sample='${sample}' type='${type}' db='${db}'"
  echo "    QIDs: $(basename "$qids")"
  echo "    OUT : $(basename "$out_fa")"

  if [[ -s "$out_fa" ]]; then
    echo "    FASTA exists, skipping build."
    return 0
  fi
  if [[ ! -s "$qids" ]]; then
    echo "    SKIP: Missing qids: $qids"
    return 1
  fi

  local src_fa
  if ! src_fa="$(find_source_fasta "$type" "$sample")"; then
    echo "    SKIP: No source FASTA found for sample '${sample}' in ${type} dir (.fasta/.fa/.fna)."
    return 1
  fi
  echo "    SRC : ${src_fa}"

  # Clean qids: remove CR, keep first token, strip [], trim trailing spaces, unique
  local tmp_qids="${qids}.clean"
  tr -d $'\r' < "$qids" \
    | awk '{print $1}' \
    | sed 's/\[\]//' \
    | sed 's/[[:space:]]\+$//' \
    | awk 'length($0)>0' \
    | sort -u > "$tmp_qids"

  # Normalize FASTA headers robustly (first token only; strip [])
  local tmp_norm="${src_fa}.norm"
  if ! seqkit fx2tab -nl "$src_fa" \
    | awk -F'\t' '{print $1"\t"$2}' \
    | awk -F'\t' '{id=$1; sub(/ .*/,"",id); gsub(/\[\]/,"",id); print id"\t"$2}' \
    | seqkit tab2fx > "$tmp_norm"; then
    echo "    ERROR: robust normalization failed for '${src_fa}'" >&2
    rm -f "$tmp_qids"
    return 1
  fi

  # Extract with seqkit grep
  local tmp="${out_fa}.tmp"
  if ! seqkit grep -f "$tmp_qids" "$tmp_norm" \
      | seqkit seq -w 0 > "$tmp"; then
    echo "    ERROR: seqkit grep failed for normalized '${tmp_norm}' with '${tmp_qids}'" >&2
    rm -f "$tmp" "$tmp_qids" "$tmp_norm"
    return 1
  fi
  rm -f "$tmp_qids" "$tmp_norm"

  if [[ ! -s "$tmp" ]]; then
    echo "    SKIP: No sequences extracted (post-normalization)."
    rm -f "$tmp"
    return 1
  fi

  mv -f "$tmp" "$out_fa"
  echo "    Built $(basename "$out_fa") with $(seqkit stats "$out_fa" | tail -n 1 | awk '{print $5}') sequences"
  return 0
}

# 1) Build per-TSV FASTAs
made_any_fasta=0
for blast_file in "${tsvs[@]}"; do
  base_name=$(basename "$blast_file" .tsv)
  if build_fasta_for_base "$base_name"; then
    [[ -s "${HCH_DIR}/${base_name}.fna" ]] && made_any_fasta=1
  fi
done

# 2) Assemble all built .fna files into ASSEMBLED_DIR
echo "= Assembling .fna results into output directory: $ASSEMBLED_DIR ="
for fna in ${HCH_DIR}/*.fna; do
  if [[ -s "$fna" ]]; then
    cp -f "$fna" "$ASSEMBLED_DIR/"
    echo "  Copied: $(basename "$fna") to $ASSEMBLED_DIR"
  fi
done

echo "= Phase 2: Coverage analysis ="

# Ensure coverage_analysis.py exists
if [[ ! -f "coverage_analysis.py" ]]; then
  echo "ERROR: coverage_analysis.py not found in current directory." >&2
  exit 1
fi

# 3) Run coverage (write directly under MASTER_COV_DIR per-sample subdir; no nested duplication)
ran_any=0
for blast_file in "${tsvs[@]}"; do
  base_name=$(basename "$blast_file" .tsv)
  fasta_file="${HCH_DIR}/${base_name}.fna"
  out_dir="${MASTER_COV_DIR}/coverage_${base_name}"

  if [[ ! -s "$fasta_file" ]]; then
    echo "SKIP: Missing FASTA for ${base_name}: ${fasta_file}"
    continue
  fi

  mkdir -p "$out_dir"
  # Files will be written directly into out_dir:
  #   ${base_name}_coverage_stats.csv
  #   ${base_name}_detailed_alignments.csv
  #   ${base_name}_coverage_analysis.png
  if [[ -s "${out_dir}/${base_name}_coverage_stats.csv" && -s "${out_dir}/${base_name}_detailed_alignments.csv" ]]; then
    echo "Coverage already present for ${base_name}, skipping."
    continue
  fi

  echo "Running coverage_analysis.py for ${base_name}"
  if ! python coverage_analysis.py --blast "$blast_file" --fasta "$fasta_file" --output_dir "$out_dir"; then
    echo "ERROR: coverage_analysis.py failed for ${base_name}" >&2
    continue
  fi
  ran_any=1
done

if [[ "$ran_any" -eq 0 ]]; then
  echo "WARNING: Coverage not run for any dataset (no matching FASTAs built or all present)."
fi

# 4) Cleanup legacy nested folders inside MASTER_COV_DIR (flatten)
echo "= Cleanup: flatten nested coverage directories (if any) ="
shopt -s nullglob
for d in "${MASTER_COV_DIR}"/coverage_*; do
  [[ -d "$d" ]] || continue
  base="${d#${MASTER_COV_DIR}/coverage_}"
  inner="${d}/${base}"
  if [[ -d "$inner" ]]; then
    echo "Flattening $d"
    find "$inner" -maxdepth 1 -type f -name "${base}_*" -exec mv -f {} "$d/" \;
    rmdir "$inner" 2>/dev/null || true
  fi
done

# 5) Success checkpoint
if ls "${MASTER_COV_DIR}"/coverage_*/*_coverage_stats.csv >/dev/null 2>&1 || [[ "$made_any_fasta" -eq 1 || "$ran_any" -eq 1 ]]; then
  echo "=== Phase 2 complete (outputs present under ${MASTER_COV_DIR}) ==="
  mkdir -p .checkpoints
  : > .checkpoints/phase2_coverage.done
else
  echo "=== Phase 2 completed with no outputs; not writing checkpoint ==="
fi
