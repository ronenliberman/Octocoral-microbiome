# Where the data came from

This file records who produced each dataset in the repository and where the
public copy lives, so that the data availability statement and any reuse can be
checked against a record rather than memory. `RULES.md` covers how each dataset
is processed.

**The repository's [MIT License](LICENSE) covers the code.** The sequence data
are public at NCBI under the accessions below and carry NCBI's terms.

Nothing below is inferred. A cell says what a document in this repository or the
published paper says, or it says that the fact is still to be recorded.

## The datasets

| Files | Contributor | Public copy | Notes |
|---|---|---|---|
| `raw files/table-NoEMC-PEE2.qza`, `taxonomy-PEE2.qza`, `dna-sequences-PEE2.fasta` | Ronen Liberman (sampling and every QIIME2/DADA2 step) | Raw reads: NCBI BioProject **PRJNA1373357** | 16S rRNA V4 (515F/806R), Illumina MiSeq 2 × 300 bp, three sampling years merged. Sequenced in house in the Lopez Lab, Nova Southeastern University. Taxonomy: SILVA naive Bayes classifier. |
| `raw files/metadata_all_edited.csv` | Ronen Liberman | with the BioProject samples | One row per sample: host (`SPP`), site, year. Site `DC2` is recoded to `DC3` in the scripts, not in the file. |
| `blast/` | NCBI BLAST and Entrez searches run by Ronen Liberman, March–October 2026 | the NCBI databases at those dates | Dates, databases and request IDs are in `blast/README.md`. |
| `phylogeny/bacteria/` | Ronen Liberman (round 3; MAFFT v7.526, IQ-TREE v3.1.4) | reference sequences are NCBI RefSeq records listed in `script/Editor_round3/00_setup.R` | Built from `raw files/` and the references. |
| Host mtMutS sequences (Fig S2; not in this repository) | Matías Gómez-Corrales | GenBank: *B. asbestinum* PX693362–PX693378; *E. caribaeorum* PX704549–PX704570; *E. flexuosa* PX704571–PX704583; *M. muricata* PX712868–PX712887 | Trees drawn by `phylogeny/NCRI_Octocoral_Host_Phylogenies_v2.Rmd`. |
| `maps/` | Matías Gómez-Corrales | | Site map (Fig 1). Site coordinates are in Table S1 of the paper. |

## Collection permits

Octocorals were collected under a Florida Fish and Wildlife Conservation
Commission Fishing License and a State of Florida Special Activity License. The
Methods give the licence as 24-2666-SRP; the Acknowledgements give
SAL-25-2666-SRP. The paper's Methods, "Sample collection", is the reference.
