
![WhatsApp Image 2025-08-21 at 14 59 16](https://github.com/user-attachments/assets/6e5e14fd-ab8f-4ed9-9e0b-24db1f54eb31)




       **viRNAtrap_project (non-canonical antigen prediction)**
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

                  
                  **Automated Pipeline**
       ┌─────────────────────────────────────────────────────────────┐
       │                    PHASE 1                                  │
       │              Extract 90% Identity Hits                      │
       │  • Filter BLAST results (≥80% identity, ≥50bp alignment,    |
       |query_len >= 20 && pident >= 80 && AL >= 50 && mismatch <= 5 | 
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

<img width="4200" height="3000" alt="microbial_coverage_plot_light" src="https://github.com/user-attachments/assets/551c3056-8e9f-4d1b-9ff1-1c2b9db83f17" />

<img width="3000" height="2400" alt="4_category_stacked_bar" src="https://github.com/user-attachments/assets/91e9af8a-104e-4e19-8bab-e39ef69ed727" />

<img width="4800" height="3600" alt="10_comprehensive_overview" src="https://github.com/user-attachments/assets/db00d8c6-5219-4389-a83d-a12ee3053725" />

<img width="3600" height="2100" alt="violin_coverage_by_group" src="https://github.com/user-attachments/assets/35e9922d-5262-4bc9-b069-acccdb914db2" />

<img width="7200" height="4800" alt="comprehensive_bacteria_cancer_analysis_classified" src="https://github.com/user-attachments/assets/e1b6671a-26e9-42c7-ab52-1b7b19e79816" />
            
<img width="4800" height="3000" alt="comprehensive_virus_cancer_analysis_classified" src="https://github.com/user-attachments/assets/141a5ab7-b5eb-4917-8137-076da74acd9f" />

<img width="4800" height="3000" alt="comprehensive_parasite_cancer_analysis_classified" src="https://github.com/user-attachments/assets/78c5169f-7114-4271-940b-d7f8d2064095" />

<img width="4800" height="3000" alt="comprehensive_parasite_cancer_analysis_unclassified" src="https://github.com/user-attachments/assets/2359e5de-c937-4ef9-ade2-f9e8c068aa19" />

<img width="4800" height="3000" alt="comprehensive_fungi_cancer_analysis_classified" src="https://github.com/user-attachments/assets/63e77459-6219-40a4-8b66-27394d85f6f8" />

<img width="4800" height="3000" alt="comprehensive_fungi_cancer_analysis_unclassified" src="https://github.com/user-attachments/assets/cf6b2d15-a57f-41c3-95cc-b24cf3485315" />


