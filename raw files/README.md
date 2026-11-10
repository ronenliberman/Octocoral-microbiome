# raw files

QIIME2 outputs and sample metadata that every analysis script starts from.
**Read-only** (AGENTS.md §2): never edit or re-save these files. This folder
holds no processed tables and no figures; the scripts write those to `output/`.

| File | What it is |
|---|---|
| `table-NoEMC-PEE2.qza` | DADA2 feature table (ASV counts per sample), chloroplast and mitochondria removed |
| `taxonomy-PEE2.qza` | SILVA naive Bayes taxonomy for each ASV, with classifier confidence |
| `dna-sequences-PEE2.fasta` | Representative sequence of each ASV (exported from `rep-seqs-PEE2.qza`) |
| `metadata_all_edited.csv` | Sample metadata: host (`SPP`), site, year |

"PEE2" marks the QIIME2 run that every published analysis uses. Raw reads are
in NCBI BioProject PRJNA1373357.

The earlier `table.qza`, `table-NoEMC.qza` and `taxonomy.qza` (from an older
QIIME2 run, used for the first-submission trees) were removed in November 2026;
no published result uses them. They remain in the git history.
