# viRNAtrap_project (non-canonical antigen prediction)

Raw RNA-seq Data
       ↓
[nfcore-rnaseq Processing]
       ↓
Unmapped RNA Sequences
       ↓
[Virnatrap AI Tool]
       ↓
    contigs
       ↓
   VecScreen
       ↓
Vector-free Sequences
       ↓
[Kraken2 Classification]
(fungi, viral, protozoa, bacteria databases)
       ↓
├─── Classified FASTA ────┐
│                         │
└─── Unclassified FASTA ──┘
        ↓
[BLAST Analysis against fungi, bacteria, parasites, virus databases]
        ↓
  BLAST Results TSVs

                  Automated Pipeline
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 1                                  │
│              Extract 90% Identity Hits                      │
│  • Filter BLAST results (≥80% identity, ≥50bp alignment)    │
│  • Generate high-confidence hits TSV and FASTA              │
└─────────────────────────────────────────────────────────────┘
       ↓ (dependency: afterok)
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 2                                  │
│                Coverage Analysis                            │
│  • Calculate per-query coverage statistics                  │
│  • Generate coverage plots and detailed alignments          │
└─────────────────────────────────────────────────────────────┘
       ↓ (dependency: afterok)
┌─────────────────────────────────────────────────────────────┐
│                   PHASE 2B                                  │
│                Coverage Filtering                           │
│  • Filter sequences by coverage thresholds (≥70%)           │
│  • Generate coverage-filtered FASTA files                   │
└─────────────────────────────────────────────────────────────┘
       ↓ (dependency: afterok)
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 3                                  │
│            Kraken2 ∩ BLAST Intersection                     │
│  • Find sequences classified by both tools                  │
│  • Generate joint classification results                    │
└─────────────────────────────────────────────────────────────┘
       ↓ (dependency: afterok)
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 4                                  │
│                Length Filtering (20bp-97bp & 97bp-200bp)    │
│  • Filter sequences by minimum length requirements          │
│  • Generate length-filtered FASTA files                     │
└─────────────────────────────────────────────────────────────┘
       ↓ (dependency: afterok)
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 5                                  │
│              Functional Annotation (DIAOMAND/pfam/cog)      │
│  • Perform functional analysis of filtered sequences        │
│  • Generate annotation results                              │
└─────────────────────────────────────────────────────────────┘
       ↓ (dependency: afterok)
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 6                                  │
│            Comprehensive Analysis (R)                       │
│  • Statistical analysis and visualization                   │
│  • Generate final reports and results                       │
└─────────────────────────────────────────────────────────────┘
                           ↓
                Final Results & Reports
