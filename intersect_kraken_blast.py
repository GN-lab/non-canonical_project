#!/usr/bin/env python3
# SBATCH --job-name=phase3
# SBATCH --partition=compute
# SBATCH --cpus-per-task=4
# SBATCH --mem-per-cpu=8042
# SBATCH --time=24:00:00
# SBATCH --output=logs/phase3_%j.log

import argparse
from pathlib import Path
from collections import defaultdict
import sys
import os
import tempfile
import csv
import re

def log(msg: str):
    print(msg, flush=True)

def atomic_write_text(path: Path, content_iterable):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(path.parent)) as tmp:
        tmp_name = tmp.name
        for line in content_iterable:
            tmp.write(line)
    os.replace(tmp_name, path)

def canon_id(s: str) -> str:
    """Consistent ID canonicalization - matches Phase 4 approach"""
    # Remove various special characters and take first token
    return re.sub(r'[\[\]|()]', '', str(s)).split()[0]

def sample_stem_from_base(base_stem: str) -> str:
    # Reduce a base like:
    #   ES01RB.unmapped_1_classified_bact_90pct
    #   ES01RB.unmapped_1_unclassified_fungi_90pct
    # to the Kraken key:
    #   ES01RB.unmapped_1
    m = re.match(r'^(.*?\.unmapped_\d+)', base_stem)
    return m.group(1) if m else base_stem

def extract_kraken_names(kraken_named_file: Path):
    """Extract Kraken classifications with proper ID parsing"""
    tax_map = {}
    with open(kraken_named_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            classification = parts[0]   # 'C' or 'U'
            seqid_raw = parts[1]        # FIXED: Get the second column (sequence ID)
            seqid = canon_id(seqid_raw)
            classified = (classification == 'C')
            tax_map[seqid] = (classified, parts)
    return tax_map

def parse_blast_hits(blast_file: Path):
    """Parse BLAST hits with proper column indexing"""
    blast_map = defaultdict(list)
    with open(blast_file) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            # BLAST outfmt 6: need at least qseqid (0) and sseqid (1)
            if len(row) < 2:
                continue
            qseqid_raw = row[0]    # FIXED: First column is qseqid
            sseqid = row[1]        # FIXED: Second column is sseqid
            qseqid = canon_id(qseqid_raw)
            blast_map[qseqid].append(sseqid)
    return blast_map

def extract_microbial_group(filename_stem: str):
    """Extract microbial group from filename"""
    parts = filename_stem.lower().split('_')
    for part in ('bact', 'virus', 'fungi', 'parasite'):
        if part in parts:
            return part
    return 'unknown'

def debug_id_matching(kraken_ids, blast_ids, sample_name):
    """Debug function to check ID matching issues"""
    log(f"DEBUG {sample_name}: Kraken IDs: {len(kraken_ids)}, BLAST IDs: {len(blast_ids)}")
    
    # Check first few IDs from each set
    kraken_sample = list(kraken_ids)[:3] if kraken_ids else []
    blast_sample = list(blast_ids)[:3] if blast_ids else []
    log(f"DEBUG {sample_name}: Kraken sample IDs: {kraken_sample}")
    log(f"DEBUG {sample_name}: BLAST sample IDs: {blast_sample}")
    
    # Check if any IDs match the pattern we expect
    if kraken_ids and blast_ids:
        common = kraken_ids & blast_ids
        log(f"DEBUG {sample_name}: Common IDs: {len(common)}")
        if common:
            log(f"DEBUG {sample_name}: First common ID: {list(common)[0]}")
        else:
            # Show why they might not match
            kraken_example = next(iter(kraken_ids))
            blast_example = next(iter(blast_ids))
            log(f"DEBUG {sample_name}: Kraken ID example: '{kraken_example}'")
            log(f"DEBUG {sample_name}: BLAST ID example: '{blast_example}'")
            log(f"DEBUG {sample_name}: After canon: '{canon_id(kraken_example)}' vs '{canon_id(blast_example)}'")

def main():
    p = argparse.ArgumentParser(description="Intersect Kraken named ids with grouped BLAST hits")
    p.add_argument("--blast_hits_dir", required=True, help="Directory of BLAST coverage filtered hits files (*.txt)")
    p.add_argument("--kraken_dir", required=True, help="Directory of Kraken named output files (*_kraken_output_with_name.txt)")
    p.add_argument("--output_dir", required=True, help="Output directory for intersection results")
    p.add_argument("--combine", action='store_true', help="Combine all results into one TSV file")
    p.add_argument("--debug", action='store_true', help="Enable debug output for ID matching")
    args = p.parse_args()

    blast_dir = Path(args.blast_hits_dir)
    kraken_dir = Path(args.kraken_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build Kraken index keyed by SAMPLE.unmapped_N (e.g., ES01RB.unmapped_1)
    kraken_index = {}
    for f in kraken_dir.glob('*_kraken_output_with_name.txt'):
        stem = f.stem.replace('_kraken_output_with_name', '')
        key = sample_stem_from_base(stem)
        kraken_index[key] = f
        if args.debug:
            log(f"Kraken file {f.name} -> key: {key}")

    all_records = []

    # Process each BLAST filtered hits file
    for blast_hit_file in sorted(blast_dir.glob('*_coverage_filtered_hits.txt')):
        base_stem = blast_hit_file.stem.replace('_coverage_filtered_hits', '')
        sample_key = sample_stem_from_base(base_stem)
        kraken_file = kraken_index.get(sample_key)
        
        if not kraken_file:
            log(f"No Kraken named file found for {base_stem} (sample key: {sample_key}), skipping")
            continue

        log(f"Processing {base_stem}")

        # Load and canonicalize IDs consistently
        kraken_names = extract_kraken_names(kraken_file)
        blast_map = parse_blast_hits(blast_hit_file)
        microbial_group = extract_microbial_group(base_stem)

        # Get ID sets for intersection
        blast_ids = set(blast_map.keys())
        kraken_ids = set(kraken_names.keys())
        intersect_ids = blast_ids & kraken_ids

        if args.debug:
            debug_id_matching(kraken_ids, blast_ids, base_stem)

        if not intersect_ids:
            log(f"No intersect IDs for {base_stem}")
            continue

        # Collect records
        local_records = []
        for qid in sorted(intersect_ids):
            classified, kname_parts = kraken_names[qid]
            blast_sids = blast_map[qid]
            kraken_name_field = "\t".join(kname_parts)
            local_records.append((qid, classified, kraken_name_field, ';'.join(blast_sids), microbial_group))

        if args.combine:
            all_records.extend(local_records)
        else:
            def tsv_lines():
                yield "query_id\tkraken_classified\tkraken_name\tblast_subject_ids\tmicrobial_group\n"
                for r in local_records:
                    yield "\t".join(map(str, r)) + "\n"

            out_file = output_dir / f"{base_stem}_kraken_blast_intersect.tsv"
            atomic_write_text(out_file, tsv_lines())
            log(f"Wrote {len(local_records)} intersection results for {base_stem} to {out_file}")

    # If combine mode, write single file
    if args.combine and all_records:
        def tsv_lines():
            yield "query_id\tkraken_classified\tkraken_name\tblast_subject_ids\tmicrobial_group\n"
            for r in all_records:
                yield "\t".join(map(str, r)) + "\n"

        combined_file = output_dir / "all_intersections.tsv"
        atomic_write_text(combined_file, tsv_lines())
        log(f"Wrote {len(all_records)} combined intersection results to {combined_file}")

if __name__ == "__main__":
    main()
