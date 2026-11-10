# blast

Saved NCBI searches used to verify the Endozoicomonadaceae and
Mycoplasmataceae ASVs (RULES.md, "ASV verification"). NCBI's databases change,
so these files are what makes the published verification reproducible.
**Read-only** (AGENTS.md §2). The BLAST scripts reuse them and query NCBI only
when a file is missing.

| File | Search | Written by | Read by |
|---|---|---|---|
| `Endo_blast_results.csv` | Endozoicomonadaceae ASVs vs NCBI 16S type strains, web BLAST hit table (March 2026) | NCBI website | `01` |
| `Myco_blast_results.csv` | Mycoplasmataceae ASVs vs NCBI 16S type strains, web BLAST hit table (March 2026) | NCBI website | `01` |
| `Mycoplasma_master_table.csv`, `Mycoplasma_master_BLAST_annotated.csv` | R1 *Mycoplasma* master tables with the Myco_1–10 labels and best type-strain hits (August 2026) | `script/archive/Mycoplasma_Phylogenetic_Analysis_UPDATED.R` | `02`, `01b_Myco` |
| `V8941X58016-Alignment-HitTable.csv` | Web BLAST of 58 ASVs from an earlier ANCOM-BC2 result set (March 2026); kept for the record, superseded by `ANCOMBC2_r3_16S_hits.csv` | NCBI website | none |
| `NCBI_BLAST_poorly_classified_ASVs_summary.csv` | Poorly classified ASVs vs a local NCBI 16S database (March 2026) | `script/NCBI_BLAST_unclassified_ASVs_localBLAST.R` | `script/Octocoral_Revised_Analysis_COMPLETE.R` |
| `ANCOMBC2_r3_16S_hits.csv`, `..._SEARCH_INFO.txt` | The 26 ASVs with significant ANCOM-BC2 results vs NCBI 16S type strains, local BLAST+ 2.17.0 (7 Oct 2026) | `script/Editor_round3/01d_ANCOMBC2_typestrain_BLAST_local.ps1` | `01` |
| `TierA_topup_16S_hits_*.csv` | PEE2 ASVs missing from the March searches vs 16S type strains (30 Sep 2026) | `01b` | `01` |
| `core_nt_raw.json`, `core_nt_hits_top10.csv`, `Unverified_ASVs_nt_hits.csv` | Tier B: unverified ASVs vs `core_nt` (30 Sep 2026) | `01b` | `01`, `02` |
| `Myco_core_nt_raw.json`, `Myco_core_nt_hits_top10.csv`, `Myco_core_nt_RID.txt` | Mycoplasma ASVs failing Tier A vs `core_nt` (RID BU04GWMW014, 30 Sep 2026) | `01b_Myco` | for the record |
| `genus_check_core_nt_WEB_HitTable*.csv` | Family-level ASVs vs `core_nt`, eight web downloads (30 Sep 2026) | NCBI website | `01c` |
| `genus_check_core_nt_hits_top10.csv`, `Family_level_ASVs_nt_hits.csv` | Genus check, parsed with NCBI lineages | `01c` | `01` |
| `NCBI_lineage_cache.rds` | NCBI taxonomy lineage for each hit accession | `01` | `01` |
| `BLAST_RIDs.txt` | NCBI request IDs of the remote searches | `01b`, `01c` | for the record |

Script numbers refer to `script/Editor_round3/`.
