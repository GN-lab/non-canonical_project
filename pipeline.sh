#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --output=/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/output/logs/megahit_%j.out
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8042

echo "Python location: $(which python3)"

echo "=== COMPLETE VIRAL ANALYSIS PIPELINE ==="
echo "Sample: ES01RB"
echo "Date: $(date)"

# Step 1: MEGAHIT paired-end assembly
echo "Step 1: Running MEGAHIT assembly..."
for merged in merged_fasta/*_merged.fastq; do
  sample=$(basename "$merged" _merged.fastq)
  outdir="megahit_assemblies/${sample}"

  echo "[$sample] Input: $merged"
  # overwrite behavior: remove old output dir if present
  if [[ -d "$outdir" ]]; then
    echo "[$sample] Removing existing output directory: $outdir"
    rm -rf "$outdir"
  fi

  megahit \
      -r "$merged" \
      --k-list 21,41,61,81,99 \
      --min-contig-len 500 \
      --presets meta-large \
      -o "$outdir" \
      -t 16

  echo "[$sample] ✓ MEGAHIT completed: $(grep -c '^>' "$outdir/final.contigs.fa") contigs generated"
done

# Step 2: Activate viRNAtrap environment and run viral detection
echo "Step 2: Switching to viRNAtrap environment..."
source /data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh
conda activate /data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/envs/env_virnatrap

# Prepare viRNAtrap input
mkdir -p virnatrap_input src
awk '/^>/{print "@"substr($0,2); next} {print; print "+"; print gensub(/./, "I", "g")}' \
    megahit_assembly/final.contigs.fa > virnatrap_input/ES01RB_contigs_unmapped.fastq

cp ../src/assemble_read_c.so src/

# Run viRNAtrap
echo "Running viRNAtrap on assembled contigs..."
python3 /data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/envs/env_virnatrap/bin/virnatrap-predict \
    --input virnatrap_input/ \
    --output virnatrap_output/ \
    --fastmode 1

echo "✓ viRNAtrap completed: $(grep -c '^>' virnatrap_output/ES01RB_contigs_contigs.txt) viral contigs detected"

# Step 3: Calculate normalized metrics
echo "Step 3: Calculating assembly metrics..."

python3 << 'EOF'
def calculate_metrics(fasta_file, name):
    lengths = []
    current_seq = ""
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_seq:
                        lengths.append(len(current_seq))
                    current_seq = ""
                else:
                    current_seq += line.strip()
            if current_seq:
                lengths.append(len(current_seq))
    except FileNotFoundError:
        print(f"File {fasta_file} not found")
        return
    
    if lengths:
        # Calculate N50
        lengths_sorted = sorted(lengths, reverse=True)
        total_length = sum(lengths_sorted)
        target = total_length * 0.5
        cumulative = 0
        n50 = 0
        
        for length in lengths_sorted:
            cumulative += length
            if cumulative >= target:
                n50 = length
                break
        
        total_contigs = len(lengths)
        contigs_1kb = len([l for l in lengths if l >= 1000])
        contigs_500bp = len([l for l in lengths if l >= 500])
        longest = max(lengths)
        
        print(f"\n=== {name} ===")
        print(f"Total contigs: {total_contigs}")
        print(f"Contigs ≥500bp: {contigs_500bp}")
        print(f"Contigs ≥1kb: {contigs_1kb}")
        print(f"N50: {n50}")
        print(f"Total bases: {total_length:,}")
        print(f"Longest contig: {longest}")
        
    else:
        print(f"No contigs found in {name}")

# Calculate metrics for both assemblies
calculate_metrics("megahit_assembly/final.contigs.fa", "MEGAHIT Assembly")
calculate_metrics("virnatrap_output/ES01RB_contigs_contigs.txt", "Viral Contigs")
EOF

echo ""
echo "✓ Complete pipeline finished!"
echo "✓ Next steps: VecScreen → Kraken2 → BLAST"
echo "✓ Results: megahit_assembly/final.contigs.fa (all contigs)"
echo "✓ Results: virnatrap_output/ES01RB_contigs_contigs.txt (viral contigs)"
