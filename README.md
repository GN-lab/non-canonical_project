




       **(non-canonical contigs prediction)**
                     Raw RNA-seq Data
                            ↓
              [nfcore-rnaseq Processing]
                            ↓
              Unmapped RNA Sequences
                            ↓
                     MEGA-HIT assembly
                            ↓
                     Assembled contigs
                            ↓
                        VecScreen
                            ↓
              Vector-free cleaned fasta
                              ↓
       [BLAST Analysis against fungi, bacteria, parasites, virus, human databases]
                             ↓
                     BLAST Results TSVs
                             ↓
                  Filter out the human hits
                             ↓
                   Filtered BLAST results
                             ↓
      Re-mapping fatsa sequence to original fasta using BWA

                  
                          **Automated Pipeline**
      
            
       ┌─────────────────────────────────────────────────────────────┐
       │                    PHASE 1                                  │
       │                Coverage Analysis                            │
       │  • Calculate per-query coverage statistics                  │
       └─────────────────────────────────────────────────────────────┘
                     ↓ (dependency: afterok)
       ┌─────────────────────────────────────────────────────────────┐
       │                   PHASE 1B                                  │
       │                Coverage Filtering                           │
       │  • Filter BLAST results (≥80% identity &&                   │
       │        coverage >= 70 && E-value <= 1e -5)                  │
       └─────────────────────────────────────────────────────────────┘
                     ↓ (dependency: afterok)
       ┌─────────────────────────────────────────────────────────────┐
       │                    PHASE 2                                  │
       │                Length Filtering (20bp-97bp & 97bp-200bp)    │
       │  • Filter sequences by minimum length requirements          │
       │  • Generate length-filtered FASTA files                     │
       └─────────────────────────────────────────────────────────────┘
                     ↓ (dependency: afterok)
       ┌─────────────────────────────────────────────────────────────┐
       │                    PHASE 3                                  │
       │              Functional Annotation (DIAOMAND/pfam/cog)      │
       │  • Perform functional analysis of filtered sequences        │
       │  • Generate annotation results                              │
       └─────────────────────────────────────────────────────────────┘
                     ↓ (dependency: afterok)
       ┌─────────────────────────────────────────────────────────────┐
       │                    PHASE 4                                  │
       │            Comprehensive Analysis (R)                       │
       │  • Statistical analysis and visualization                   │
       │  • Generate final reports and results                       │
       └─────────────────────────────────────────────────────────────┘
                               ↓
                     Final Results & Reports

<img width="4800" height="3600" alt="10_comprehensive_overview" src="https://github.com/user-attachments/assets/aabd266d-90ce-45ef-8354-218851d44376" />

<img width="3000" height="2400" alt="4_category_stacked_bar" src="https://github.com/user-attachments/assets/d25f1e67-b9e9-44b6-ab89-9307b2d66035" />


<img width="2400" height="1800" alt="2_sample_category_distribution" src="https://github.com/user-attachments/assets/7bc2dbbb-bc34-4745-a6fb-60354fa84a74" />

<img width="2400" height="1800" alt="1_microbial_group_distribution" src="https://github.com/user-attachments/assets/eabaf8a8-5c58-4616-a30a-3c8f36e0b2e1" />

<img width="2400" height="1800" alt="8_relative_abundance_heatmap" src="https://github.com/user-attachments/assets/00aceda7-3623-4d77-b8c4-1810ac6c7145" />

<img width="4800" height="3000" alt="comprehensive_virus_cancer_analysis" src="https://github.com/user-attachments/assets/fbd23f47-465b-4fa8-af2f-700c20ce27fa" />

<img width="7200" height="4800" alt="comprehensive_bacteria_cancer_analysis" src="https://github.com/user-attachments/assets/7b80423f-d28e-46fa-b2f7-df9cd84e7172" />

<img width="4800" height="3000" alt="comprehensive_fungi_cancer_analysis" src="https://github.com/user-attachments/assets/3b51b224-285a-4e32-b32e-e495add27705" />

<img width="4800" height="3000" alt="comprehensive_parasite_cancer_analysis" src="https://github.com/user-attachments/assets/98eaa384-769e-4143-a45a-50f27aed1b99" />
