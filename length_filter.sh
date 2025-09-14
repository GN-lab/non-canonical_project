#!/bin/bash
#SBATCH --job-name=phase4
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8042
#SBATCH --time=24:00:00
#SBATCH --output=logs/phase4_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "ERROR rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Activate the virtual environment (ensure biopython is available; pip install biopython if needed)
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$THREADS"

INPUT_DIR="coverage_filtered"
OUTPUT_DIR="length_filtered"
COV_ANALYSIS_DIR="cov_analysis"
mkdir -p "$OUTPUT_DIR" "$COV_ANALYSIS_DIR"

filter_by_length() {
    local input_fasta=$1
    local min_length=$2
    local max_length=$3
    local output_fasta=$4

    # Idempotence
    if [[ -s "$output_fasta" ]]; then
        echo "Exists: $output_fasta (skipping)"
        return
    fi

    # Inline Python with correct argv usage and atomic write
    python3 - "$input_fasta" "$min_length" "$max_length" "$output_fasta" <<'PY'
import sys, os
from Bio import SeqIO

# Expect: sys.argv[1]=input_fasta, [2]=min_len, [3]=max_len, [4]=output_fasta
inp = sys.argv[1]
min_len = int(sys.argv[2])
max_len = int(sys.argv[3])
outp = sys.argv[4]

if max_len <= 0:
    max_len = float('inf')

tmp = outp + ".tmp"
with open(tmp, 'w') as out_handle:
    for record in SeqIO.parse(inp, 'fasta'):
        L = len(record.seq)
        if min_len <= L <= max_len:
            SeqIO.write(record, out_handle, 'fasta')
os.replace(tmp, outp)
print(f"Filtered sequences written to {outp}")
PY
}

# Length filtering categories
shopt -s nullglob
inputs=( "${INPUT_DIR}"/*_coverage_filtered.fna )
if (( ${#inputs[@]} == 0 )); then
    echo "ERROR: No input *_coverage_filtered.fna files in ${INPUT_DIR}" >&2
    exit 1
else
    for fasta_file in "${inputs[@]}"; do
        base_name=$(basename "$fasta_file" _coverage_filtered.fna)
        echo "Processing $base_name..."

        # 48–97bp (adjusted min to 20 as per Phase 1 filter)
        filter_by_length "$fasta_file" 20 97 "${OUTPUT_DIR}/${base_name}_48-97bp.fna"
        # 98–200bp
        filter_by_length "$fasta_file" 98 200 "${OUTPUT_DIR}/${base_name}_98-200bp.fna"
        # >200bp (use max_length=0 to mean 'no upper bound')
        filter_by_length "$fasta_file" 201 0 "${OUTPUT_DIR}/${base_name}_long.fna"

        # Special case: fungi — copy full-length for supplementary analysis
        if [[ "$base_name" == *"fungi"* ]]; then
            out_all="${OUTPUT_DIR}/${base_name}_all_lengths.fna"
            if [[ -s "$out_all" ]]; then
                echo "Exists: $out_all (skipping copy)"
            else
                cp "$fasta_file" "$out_all"
                echo "Fungi sequences saved for comprehensive functional analysis"
            fi
        fi
    done
fi

# Extract canonical IDs for coverage summary (first token only)
echo "Extracting filtered sequence IDs for coverage summary..."
shopt -s nullglob
for fasta_file in "${OUTPUT_DIR}"/*_*.fna; do
  ids_out="${fasta_file%.fna}_ids.txt"
  if [[ -s "$ids_out" ]]; then
    echo "Exists: $ids_out (skipping)"
  else
    # Strip '>' and take only the first whitespace-delimited token (canonical ID)
    awk '/^>/ {sub(/^>/,"",$0); split($0,a," "); print a[0]}' "$fasta_file" > "${ids_out}.tmp"
    mv "${ids_out}.tmp" "$ids_out"
  fi
done

# R coverage-depth summary (replaced with Python as in original)
echo "Running python to create coverage_depth_stats.tsv..."

# Activate the virtual environment
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

python3 - <<'EOF'
import os
import glob
import subprocess
import collections
import re

def run_seqkit_cmd(cmd):
    """Run seqkit command and return output"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.stdout.strip()
    except Exception as e:
        print(f"Error running seqkit: {e}")
        return ""

input_dir = "coverage_filtered"
output_dir = "length_filtered"
cov_analysis_dir = "cov_analysis"
os.makedirs(cov_analysis_dir, exist_ok=True)

summary_file = os.path.join(cov_analysis_dir, "coverage_depth_stats.tsv")

print("=== Using SeqKit for ID extraction ===")
# Use seqkit to extract IDs from all FASTA files
for fasta_file in glob.glob(os.path.join(output_dir, "*.fna")):
    ids_file = fasta_file.replace('.fna', '_ids.txt')
    
    if os.path.getsize(fasta_file) > 0:
        # Get sequence count
        stats_cmd = f"seqkit stats {fasta_file} -T"
        stats = run_seqkit_cmd(stats_cmd)
        num_seqs = stats.split('\n')[1].split()[3] if stats else "0"
        
        print(f"Processing {os.path.basename(fasta_file)} ({num_seqs} sequences)")
        
        # Extract IDs using seqkit
        id_cmd = f"seqkit seq -n -i {fasta_file} | cut -d' ' -f1"
        ids_output = run_seqkit_cmd(id_cmd)
        
        with open(ids_file, 'w') as f:
            f.write(ids_output)
        
        num_ids = len(ids_output.split('\n')) if ids_output else 0
        print(f"  Extracted {num_ids} IDs")
    else:
        print(f"Empty: {os.path.basename(fasta_file)}")
        open(ids_file, 'w').close()

print("=== Creating coverage summary ===")
with open(summary_file, 'w') as out_f:
    out_f.write("Sample\tMicrobialGroup\tqseqid\tNumAlignments\n")
    
    for blast_file in glob.glob(os.path.join(input_dir, "*_coverage_filtered_hits.txt")):
        base_name = os.path.basename(blast_file).replace('_coverage_filtered_hits.txt', '')
        
        # Determine microbial group
        base_lower = base_name.lower()
        if 'virus' in base_lower:
            group = 'VIRUS'
        elif 'bact' in base_lower:
            group = 'BACT'
        elif 'parasite' in base_lower:
            group = 'PARASITE'
        elif 'fungi' in base_lower:
            group = 'FUNGI'
        else:
            group = 'OTHER'
        
        print(f"Processing {base_name} ({group})")
        
        # Collect all filtered IDs for this sample
        filtered_ids = set()
        pattern = os.path.join(output_dir, f"{base_name}_*_ids.txt")
        for id_file in glob.glob(pattern):
            with open(id_file) as f:
                for line in f:
                    if line.strip():
                        filtered_ids.add(line.strip())
        
        print(f"  Found {len(filtered_ids)} filtered IDs")
        if not filtered_ids:
            continue
        
        # Count alignments in BLAST file
        alignment_count = collections.defaultdict(int)
        blast_lines = 0
        
        try:
            with open(blast_file) as f:
                for line in f:
                    if line.strip() and not line.startswith(('qseqid', '#')):
                        blast_lines += 1
                        parts = line.strip().split('\t')
                        if parts:
                            qseqid = parts[0]
                            # Clean up the ID
                            clean_id = re.sub(r'[\[\]|()]', '', qseqid).split()[0]
                            
                            if clean_id in filtered_ids:
                                alignment_count[clean_id] += 1
        except Exception as e:
            print(f"  Error reading {blast_file}: {e}")
            continue
        
        print(f"  Processed {blast_lines} BLAST lines")
        
        # Write results
        for qseqid, count in alignment_count.items():
            out_f.write(f"{base_name}\t{group}\t{qseqid}\t{count}\n")
        
        print(f"  Added {len(alignment_count)} entries")

print(f"=== Summary written to {summary_file} ===")
# Show stats
if os.path.exists(summary_file):
    with open(summary_file) as f:
        lines = f.readlines()
        print(f"Total entries: {len(lines) - 1}")  # exclude header
        if len(lines) > 1:
            print("First few entries:")
            for line in lines[1:6]:  # first 5 data lines
                print(line.strip())
else:
    print("No summary file created")

print("Done!")

print("Generating per microbial group tsv with subject ids")
import os
import glob
import re
import csv

COVERAGE_DIR = "coverage_filtered"
OUT_DIR = "cov_analysis"
BLAST_DIR = "../filtered_blast_results"  # Single directory for filtered BLAST TSVs

GROUP_MAP = {
    "virus": "virus",
    "viruses": "virus",
    "bact": "bact",
    "bacteria": "bact",
    "bacterial": "bact",
    "parasite": "parasite",
    "parasites": "parasite",
    "fungi": "fungi",
    "fungus": "fungi",
}

def canon_id(s: str) -> str:
    if not s:
        return ""
    s = re.sub(r'[\[\]|()]', '', str(s))
    return s.split()[0].lower().strip()

def parse_sample(sample_basename: str):
    stem = sample_basename
    stem = re.sub(r"_coverage_filtered_hits$", "", stem)
    stem = re.sub(r"_90pct(_pass)?$", "", stem)

    # Updated parsing for new filenames (e.g., SAMPLE_vs_GROUP_filtered)
    m = re.match(
        r"^(?P<prefix>[^_]+)_vs_(?P<group>[^_]+)$",
        stem,
        re.IGNORECASE,
    )
    if not m:
        return None

    prefix = m.group("prefix")
    group_raw = m.group("group").lower()
    group_norm = GROUP_MAP.get(group_raw, group_raw)
    return prefix, group_norm

def build_blast_path(prefix: str, group: str) -> str | None:
    fname = f"{prefix}_vs_{group}_filtered.tsv"
    path = os.path.join(BLAST_DIR, fname)
    return path if os.path.isfile(path) else None

def read_keep_qids(coverage_file: str) -> set[str]:
    keep = set()
    with open(coverage_file) as f:
        reader = csv.reader(f, delimiter="\t")
        first = True
        for row in reader:
            if not row:
                continue
            if first:
                first = False
                if row[0].strip().lower() == "qseqid":
                    continue
            qid = canon_id(row[0])
            if qid:
                keep.add(qid)
    return keep

def stream_blast_rows(blast_path: str, keep_qids: set[str]):
    """
    Stream BLAST rows, only returning those with qids in keep_qids
    """
    # CRITICAL: If keep_qids is empty, don't yield anything
    if not keep_qids:
        return
    
    with open(blast_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 12:
                continue
            qid = canon_id(row[0])
            if qid in keep_qids:
                yield qid, row[1], row[2], row[10], row[11], row[3]

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Prepare output files for each group (no classified/unclassified distinction)
    groups = ["FUNGI", "BACT", "VIRUS", "PARASITE"]
    writers = {}
    handles = {}
    
    for grp in groups:
        path = os.path.join(OUT_DIR, f"{grp}_coverage_hits.tsv")
        fh = open(path, "w", newline="")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Sample", "QueryID", "SubjectID", "Identity", "Evalue", "Bitscore", "AlignmentLength"])
        writers[grp] = w
        handles[grp] = fh

    total = 0
    processed_files = 0
    skipped_files = 0
    
    try:
        for coverage_path in sorted(glob.glob(os.path.join(COVERAGE_DIR, "*_coverage_filtered_hits.txt"))):
            sample_base = os.path.basename(coverage_path).removesuffix(".txt")
            parts = parse_sample(sample_base)
            if not parts:
                skipped_files += 1
                continue
            
            prefix, group_norm = parts
            out_grp = group_norm.upper()
            
            if out_grp not in writers:
                print(f"Skipping {sample_base} - group {out_grp} not in output groups")
                skipped_files += 1
                continue

            print(f"Processing {sample_base}")
            
            keep_qids = read_keep_qids(coverage_path)
            if not keep_qids:
                print(f"  No qids to keep, skipping")
                skipped_files += 1
                continue

            blast_path = build_blast_path(prefix, group_norm)
            if not blast_path:
                print(f"  BLAST file not found")
                skipped_files += 1
                continue

            matched = 0
            # This is where the filtering happens - only rows with qids in keep_qids are processed
            for qid, sseqid, pident, evalue, bitscore, length in stream_blast_rows(blast_path, keep_qids):
                writers[out_grp].writerow([sample_base, qid, sseqid, pident, evalue, bitscore, length])
                matched += 1
                total += 1

            print(f"  Wrote {matched} rows")
            processed_files += 1

    finally:
        for fh in handles.values():
            fh.close()

    print(f"\n=== Summary ===")
    print(f"Processed {processed_files} files with data")
    print(f"Skipped {skipped_files} files")
    print(f"Wrote {total} total filtered rows to group TSVs in {OUT_DIR}")

if __name__ == "__main__":
    main()
EOF
