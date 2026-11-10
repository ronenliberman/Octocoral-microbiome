# Where the data came from

This file records who produced each dataset in the repository and where the
public copy lives, so that the data availability statement and any reuse can be
checked against a record rather than memory. `RULES.md` covers how each dataset
is processed.

**The repository's [MIT License](LICENSE) covers the code.** The sequence data
are public at NCBI under the accessions below and carry NCBI's terms.

Nothing below is inferred. A cell says what a document in this repository or the
published paper says, or it says that the fact is still to be recorded. Fill an
empty cell from a record (a facility report, a permit, an email) and name the
source.

## The datasets

| Files | Contributor | Public copy | Notes |
|---|---|---|---|
| `raw files/table-NoEMC-PEE2.qza`, `taxonomy-PEE2.qza`, `dna-sequences-PEE2.fasta` | Ronen Liberman (sampling); who ran each step from extraction to QIIME2/DADA2: to record | Raw reads: NCBI BioProject **PRJNA1373357** | 16S rRNA V4 (515F/806R), Illumina MiSeq 2 × 300 bp, three sampling years merged. Taxonomy: SILVA naive Bayes classifier. Sequencing facility: to record. |
| `raw files/metadata_all_edited.csv` | Ronen Liberman | with the BioProject samples | One row per sample: host (`SPP`), site, year. Site `DC2` is recoded to `DC3` in the scripts, not in the file. |
| `blast/` | NCBI BLAST and Entrez searches run by Ronen Liberman, March–October 2026 | the NCBI databases at those dates | Dates, databases and request IDs are in `blast/README.md`. |
| `phylogeny/bacteria/` | Ronen Liberman (round 3; MAFFT v7.526, IQ-TREE v3.1.4) | reference sequences are NCBI RefSeq records listed in `script/Editor_round3/00_setup.R` | Built from `raw files/` and the references. |
| Host mtMutS sequences (Fig S2; not in this repository) | Matías Gómez-Corrales | GenBank: *B. asbestinum* PX693362–PX693378; *E. caribaeorum* PX704549–PX704570; *E. flexuosa* PX704571–PX704583; *M. muricata* PX712868–PX712887 | Trees drawn by `phylogeny/NCRI_Octocoral_Host_Phylogenies_v2.Rmd`. |
| `maps/` | Matías Gómez-Corrales | | Site map (Fig 1). Site coordinates are in Table S1 of the paper. |

## Still to record

- The sequencing facility and run dates for the three MiSeq runs.
- The collection permit number(s).
