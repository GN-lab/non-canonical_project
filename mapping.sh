#!/usr/bin/env bash
#SBATCH --job-name=BWA_mapping
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8042
#SBATCH --time=48:00:00
#SBATCH --output=logs/BWA_%j.log

set -euo pipefail

module load BWA/0.7.17
module load SAMtools/1.11

# Define directories
UNMAPPED_DIR_ESO="../../nextflow_rnaseq/output_ESO/star_salmon/unmapped"
UNMAPPED_DIR_HT="../../nextflow_rnaseq/HTresults/star_salmon/unmapped"
UNMAPPED_DIR_hartwig="../../nextflow_rnaseq/output_hartwig/star_salmon/unmapped"
CLEANED_FASTA_DIR="../output/cleaned_fastas"
OUTPUT_DIR="read_mapping_results"

mkdir -p "$OUTPUT_DIR" logs

echo "=== Starting BWA mapping for all samples ==="

# Find all cleaned fasta files
for ASSEMBLY in "$CLEANED_FASTA_DIR"/*.clean.fasta; do
    # Extract sample name from assembly filename
    SAMPLE=$(basename "$ASSEMBLY" .clean.fasta)
    
    echo "Processing sample: $SAMPLE"
    
    # Determine which unmapped directory to use based on sample prefix
    if [[ "$SAMPLE" == ES* ]]; then
        UNMAPPED_DIR="$UNMAPPED_DIR_ESO"
        echo "Using ESO directory for $SAMPLE"
    elif [[ "$SAMPLE" == HT* ]]; then
        UNMAPPED_DIR="$UNMAPPED_DIR_HT"
        echo "Using HT directory for $SAMPLE"
    else
        UNMAPPED_DIR="$UNMAPPED_DIR_hartwig"
        echo "Using hartwig directory for $SAMPLE"
    fi
    
    # Define read files
    R1_READS="$UNMAPPED_DIR/${SAMPLE}.unmapped_1.fastq"
    R2_READS="$UNMAPPED_DIR/${SAMPLE}.unmapped_2.fastq"
    
    # Check if read files exist
    if [[ ! -f "$R1_READS" || ! -f "$R2_READS" ]]; then
        echo "Warning: Read files not found for $SAMPLE"
        echo "  Expected R1: $R1_READS"
        echo "  Expected R2: $R2_READS"
        echo "  Skipping..."
        continue
    fi
    
    echo "Found read files:"
    echo "  R1: $R1_READS"
    echo "  R2: $R2_READS"
    
    echo "=== Mapping $SAMPLE reads back to assembly ==="
    
    # 1. Index the assembly if needed
    if [[ ! -f "${ASSEMBLY}.bwt" ]]; then
        echo "Step 1: Indexing assembly for $SAMPLE..."
        bwa index "$ASSEMBLY"
    else
        echo "Step 1: Assembly already indexed for $SAMPLE, skipping."
    fi
    
    # 2. Map paired-end reads and sort
    echo "Step 2: Mapping paired reads for $SAMPLE..."
    bwa mem -t 16 "$ASSEMBLY" "$R1_READS" "$R2_READS" \
        | samtools sort -@ 4 -o "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam"
    
    # 3. Index BAM for analysis
    echo "Step 3: Indexing BAM for $SAMPLE..."
    samtools index "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam"
    
    # 4. Generate mapping statistics
    echo "Step 4: Generating statistics for $SAMPLE..."
    
    samtools flagstat "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam" \
        > "$OUTPUT_DIR/${SAMPLE}.flagstat.txt"
    
    samtools idxstats "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam" \
        > "$OUTPUT_DIR/${SAMPLE}.idxstats.tsv"
    
    samtools depth "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam" \
        > "$OUTPUT_DIR/${SAMPLE}.depth.tsv"
    
    # Quick summary
    echo "=== Quick Summary for $SAMPLE ==="
    total_reads=$(samtools view -c "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam")
    mapped_reads=$(samtools view -c -F 4 "$OUTPUT_DIR/${SAMPLE}_reads_vs_contigs.bam")
    mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
    
    echo "Total alignments: $total_reads"
    echo "Mapped alignments: $mapped_reads"
    echo "Mapping rate:      $mapping_rate%"
    
    # Contigs with read support
    contigs_with_reads=$(awk '$3 > 0 {count++} END {print count+0}' "$OUTPUT_DIR/${SAMPLE}.idxstats.tsv")
    total_contigs=$(awk 'END {print NR}' "$OUTPUT_DIR/${SAMPLE}.idxstats.tsv")
    echo "Contigs with support: $contigs_with_reads / $total_contigs"
    echo ""
    
done

echo "=== All samples processed ==="
