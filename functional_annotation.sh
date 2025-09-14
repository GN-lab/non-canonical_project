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

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$THREADS"
echo "Using THREADS=$THREADS"

BLAST_BIN="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Blast/bin"
export PATH="$BLAST_BIN:$PATH"

# Activate the virtual environment with bio installed
source "/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/Neojuction_pred/venv/bin/activate"

INPUT_DIR="length_filtered"
OUTPUT_DIR="functional_analysis"
mkdir -p "$OUTPUT_DIR"

# Database paths - set to the actual path for your system
DIAMOND_NR="/home/csalas/csalas_rds/gaurav_rds/virnatrap/has_output/ref_dbs/databases/nr.dmnd"
PFAM_DB="/home/csalas/csalas_rds/gaurav_rds/virnatrap/has_output/ref_dbs/databases/Pfam-A.hmm"
COG_DB="/home/csalas/csalas_rds/gaurav_rds/virnatrap/has_output/ref_dbs/databases/cog.dmnd"

annotate_sequences() {
    local input_fasta=$1
    local sample_name
    sample_name=$(basename "$input_fasta" .fna)
    local output_prefix="${OUTPUT_DIR}/${sample_name}"

    echo "Annotating $sample_name..."

    # 1) DIAMOND BLASTX against nr
    if [[ -f "$DIAMOND_NR" ]]; then
        local out_nr="${output_prefix}_nr.tsv"
        if [[ -s "$out_nr" ]]; then
            echo "  DIAMOND BLASTX already done, skipping."
        else
            echo "  Running DIAMOND BLASTX with --threads $THREADS..."
            local tmp="${out_nr}.tmp"
            diamond blastx \
                --db "$DIAMOND_NR" \
                --query "$input_fasta" \
                --out "$tmp" \
                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
                --evalue 1e-3 \
                --max-target-seqs 5 \
                --threads "$THREADS"
            mv -f "$tmp" "$out_nr"
        fi
    else
        echo "  WARNING: DIAMOND NR database not found at $DIAMOND_NR; skipping BLASTX."
    fi

    # 2) Translate sequences for protein domain analysis
    local translated="${output_prefix}_translated.faa"
    if [[ -s "$translated" ]]; then
        echo "  Translated FASTA exists, skipping translation."
    else
        echo "  Translating sequences..."
        INPUT_FASTA="$input_fasta" OUTPUT_FAA="$translated" python3 - <<'PY'

from Bio import SeqIO
import os

input_file = os.environ['INPUT_FASTA']
output_file = os.environ['OUTPUT_FAA']

def translate_all_frames(seq):
    translations = []
    # forward frames
    for frame in range(3):
        subseq = seq[frame:]
        if len(subseq) >= 3:
            prot = subseq.translate(to_stop=False)
            translations.append(str(prot))
    # reverse complement frames
    rev_comp = seq.reverse_complement()
    for frame in range(3):
        subseq = rev_comp[frame:]
        if len(subseq) >= 3:
            prot = subseq.translate(to_stop=False)
            translations.append(str(prot))
    return translations

tmp = output_file + ".tmp"
with open(tmp, 'w') as out_handle:
    for record in SeqIO.parse(input_file, 'fasta'):
        translations = translate_all_frames(record.seq)
        for j, translation in enumerate(translations):
            if len(translation.replace('*','')) >= 10:
                out_handle.write(f">{record.id}_frame_{j}\n{translation}\n")
os.replace(tmp, output_file)
PY
    fi

    # 3) HMMER search against Pfam
    if [[ -f "$PFAM_DB" && -s "$translated" ]]; then
        local pfam_tbl="${output_prefix}_pfam.tbl"
        local pfam_domtbl="${output_prefix}_pfam_domains.tbl"
        local pfam_out="${output_prefix}_pfam.out"
        if [[ -s "$pfam_tbl" && -s "$pfam_domtbl" ]]; then
            echo "  HMMER Pfam search already done, skipping."
        else
            echo "  Running HMMER against Pfam with --cpu $THREADS..."
            local tmp_tbl="${pfam_tbl}.tmp"
            local tmp_dom="${pfam_domtbl}.tmp"
            hmmsearch \
                --tblout "$tmp_tbl" \
                --domtblout "$tmp_dom" \
                -E 1e-5 \
                --cpu "$THREADS" \
                "$PFAM_DB" \
                "$translated" > "$pfam_out"
            mv -f "$tmp_tbl" "$pfam_tbl"
            mv -f "$tmp_dom" "$pfam_domtbl"
        fi
    else
        echo "  WARNING: Missing Pfam DB or translated FASTA; skipping HMMER."
    fi

    # 4) COG annotation (if configured)
    if [[ -f "$COG_DB" && -s "$translated" ]]; then
        local cog_out="${output_prefix}_cog.tsv"
        if [[ -s "$cog_out" ]]; then
            echo "  COG DIAMOND already done, skipping."
        else
            echo "  Running DIAMOND BLASTP to COG with --threads $THREADS..."
            local tmp="${cog_out}.tmp"
            diamond blastp \
                --db "$COG_DB" \
                --query "$translated" \
                --out "$tmp" \
                --outfmt 6 \
                --evalue 1e-5 \
                --max-target-seqs 1 \
                --threads "$THREADS"
            mv -f "$tmp" "$cog_out"
        fi
    fi
}


# Process all sequence files
shopt -s nullglob
for fasta_file in "$INPUT_DIR"/*.fna; do
    if [[ -s "$fasta_file" ]]; then
        annotate_sequences "$fasta_file"
    fi
done
