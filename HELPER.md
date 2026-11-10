# HELPER.md: how this project is put together

This file describes the structure: what lives where, what runs in what order,
and what each script reads and writes. **The working rules are in
[AGENTS.md](AGENTS.md)** and apply to everyone, including any assistant. **How
each dataset is handled, and why, is in [RULES.md](RULES.md).** The science and
the review history are in
[documentation/PROJECT_CONTEXT.md](documentation/PROJECT_CONTEXT.md). Update
this file in the same commit as any change to the layout.

## What the project is

16S rRNA (V4) amplicon analysis of the bacterial communities of four Florida
octocorals, the branching *Eunicea flexuosa* and *Muricea muricata* and the
encrusting *Briareum asbestinum* and *Erythropodium caribaeorum*, sampled at four
reef sites (Broward and Miami, shallow and deep) in 2023, 2024 and 2025, with
seawater controls. It is the code for Liberman, Gómez-Corrales and Lopez, *FEMS
Microbiology Ecology* (fiag124).

The work went through three rounds of review. The round-2 scripts in `script/`
do the community analysis (diversity, PERMANOVA, core microbiome, ANCOM-BC2).
The round-3 scripts in `script/Editor_round3/` verify the *Endozoicomonas* and
*Mycoplasma* ASVs against NCBI, give every ASV one stable label, rebuild the two
bacterial trees and redraw Figures 3–6 and the supplementary workbook.

## Layout

```
raw files/               QIIME2 outputs (PEE2) + sample metadata. Read-only (AGENTS.md §2).
blast/                   Saved NCBI BLAST results and taxonomy cache. Read-only (AGENTS.md §2).
phylogeny/
  bacteria/              Round-3 Endozoicomonadaceae / Mycoplasmataceae tree inputs, MAFFT
                         alignments and IQ-TREE results (Figs 5-6).
  *.Rmd, README          Host mtMutS phylogenies (Fig S2), by Matías Gómez-Corrales.
maps/                    Site map (Fig 1), by Matías Gómez-Corrales.
script/                  Round-2 analysis scripts.
  Editor_round3/         Round-3 scripts 00-17 and two PowerShell runners.
  archive/               Superseded scripts, kept for the record. Not part of the pipeline.
documentation/           PROJECT_CONTEXT.md and the round-1/2 analysis notes.
output/                  Everything the scripts write. Ignored by git; created on the first run.
.github/                 Issue and pull request templates, the attribution check workflow.
.githooks/pre-push       Local copy of the attribution check (git config core.hooksPath .githooks).
tools/check-attribution.sh  The single definition of the no-attribution rule.
```

## Paths

Every path is relative to the repository root. `here::here()` finds the root by
the `.Rproj` file, so open `Octocoral-microbiome.Rproj` (or start R anywhere in
the repository) and source a script; nothing needs editing.
`script/Editor_round3/00_setup.R` defines the folders the round-3 scripts use:

| Variable | Folder |
|---|---|
| `raw_dir` | `raw files/` |
| `blast_dir` | `blast/` |
| `mac_dir` | `phylogeny/bacteria/` (historical name: the first trees were run on a Mac) |
| `old_out_dir` | `output/` (round-2 script outputs) |
| `out_dir` | `output/Editor_round3/` |
| `supp_dir` | `output/Supplementary/` |

## Run order

### Round 2 (`script/`)

Run in one R session, in this order. The later scripts use objects that the
first one leaves in memory (`d`, `mt_rarefied`, `mt_TSS`, `mt_CLR`).

| # | Script | Reads | Writes (in `output/`) | Paper |
|---|---|---|---|---|
| 1 | `Octocoral_Revised_Analysis_COMPLETE.R` | `raw files/`, `blast/NCBI_BLAST_poorly_classified_ASVs_summary.csv` | alpha diversity, PCoA/NMDS, PERMANOVA, global ANCOM-BC2 | Fig 2B–D, Table 2 |
| 2 | `Community compsition ASV_top 20 percent per sample_MARCH17.R` | objects from 1 | stacked bar chart | Fig 2A |
| 3 | `Core_Microbiome_Revised_AllThresholds.R` | `raw files/` (standalone) | `Core_Microbiome_AllThresholds_Summary_TSS.csv` | input to round 3 |
| 4 | `Supplementary_Materials.R` | objects from 1; output of 3 and of round-3 `06c` | `output/Supplementary/` | Table S2, species PCoA panels, and the sources of the alpha-diversity, PERMANOVA and ANCOM-BC2 supplementary tables |
| – | `NCBI_BLAST_unclassified_ASVs_localBLAST.R` | `raw files/`; a local NCBI 16S database | the BLAST summary saved in `blast/` | only to redo that search |
| – | `Figure4_species_ANCOMBC_REVISION2.R` | not run directly | | run through round-3 `06c` and `06b` |

Figure 2 was assembled from the panels of scripts 1 and 2 in GIMP; round-3
script `17` then adds the panel A y-axis title.

### Round 3 (`script/Editor_round3/`)

Every script sources `00_setup.R` (paths, packages, species palette,
thresholds, reference lists, tree-plotting helpers). Run each in a fresh
session, in this order.

| # | Script | Reads | Writes | Notes |
|---|---|---|---|---|
| 1 | `06c_rerun_ANCOMBC2_PEE2.R` | `raw files/`, the model part of `Figure4_species_ANCOMBC_REVISION2.R` | `output/Editor_round3/ANCOMBC2_PEE2/` | ~5 min. Year and site effects within each host |
| 2 | `01d_ANCOMBC2_typestrain_BLAST.R` | output of 06c | `ANCOMBC2_r3_ASVs.fasta`; `blast/ANCOMBC2_r3_16S_hits.csv` | **network**. The saved hits came from the local runner `01d_..._local.ps1`. Skip to reuse `blast/` |
| 3 | `01_BLAST_validation.R` | `blast/` hit tables, lineage cache | `ASV_validation_status.csv`, Appendix 1, `Unverified_ASVs_for_nt_BLAST.fasta` | Tier A + Tier B rule (RULES.md) |
| 4 | `01b_core_nt_BLAST.R`, `01b_core_nt_BLAST_Myco.R`, `01c_core_nt_genus_upgrade.R` | outputs of 01; `blast/` | `blast/` (cached; re-query only if the cache is absent) | **network**. Then rerun 01. Skip to reuse `blast/` |
| 5 | `02_stable_ASV_IDs.R` | 01, round-2 core table, 06c, `blast/` | `ASV_stable_IDs.csv` | **network** (one NCBI taxonomy lookup). Stops if a label is duplicated |
| 6 | `03_Endo_tree_revised.R`, `04_Myco_tree_revised.R` | 02, `phylogeny/bacteria/` | `*_all_ASVs_with_status.csv`, alignment stats, `phylogeny/bacteria/*_for_tree.fasta`, `*_reference_names.csv` | Phase 1 writes the tree FASTA; phase 2 reads the alignment and tree |
| – | `11_run_alignment_and_tree.ps1` | `phylogeny/bacteria/*_for_tree.fasta` | MAFFT alignment and IQ-TREE files in `phylogeny/bacteria/` | only to rebuild the trees (MAFFT v7.526, IQ-TREE v3.1.4) |
| 7 | `05_Figure3_core_revised.R` | `raw files/`, 02 | Figure 3, `Appendix2_core_microbiome_all_thresholds.csv` | |
| 8 | `06b_run_Figure4_relabelled.R` (+ `06_Figure4_relabel.R`) | 06c, 02, the plotting part of `Figure4_species_ANCOMBC_REVISION2.R` | Figure 4 (`Figure4_ANCOMBC2_bubble_ALTERNATIVE_revision2.png`) | |
| 9 | `07_TableS8_revised.R` | `phylogeny/bacteria/` alignments, 03/04 | `TablesS8_S9_revised.xlsx` | pairwise identity, validation table |
| 10 | `08_supplement_export.R` | 01, 02, 05, 07, `raw files/` | `Supplementary_Data_Liberman_et_al.xlsx` | core bootstrap (1,000 subsamples, seed 123) |
| 11 | `16_tree_figures_v3.R`, `16b_Fig5_v3_legend_positions.R` | 03/04 outputs, `phylogeny/bacteria/` | `Figure6_Mycoplasmataceae_tree_v3.png`, `Figure5_Endozoicomonadaceae_tree_v3b.png` | final Figs 5 and 6 |
| 12 | `17_Fig2_y_axis_title.R` | `output/Editor_round3/Figure2_combined_final_GIMP.png` (place it there) | `Figure2_combined_final_ytitle.png` | final Fig 2 |

Figure 1 is `maps/Florida_map_Liberman_2025_v2.Rmd`. The host phylogenies (Fig
S2) are `phylogeny/NCRI_Octocoral_Host_Phylogenies_v2.Rmd`; their mtMutS
alignments are not in this repository (sequences are in GenBank, see
`DATA_PROVENANCE.md`).

The supplementary workbook names its sheets by the round-3 numbering (Tables
S8–S10). The published supplement renumbered them; the workbook README sheet
and the supplement captions give the final numbers.

## What is and is not tracked

Tracked: `raw files/`, `blast/`, `phylogeny/` (including the round-3 alignments
and trees), `maps/`, every script, the documentation.

Not tracked: `output/` (all figures, tables and logs the scripts write), R
session files, Office lock files, and any image file (`*.png`, `*.pdf`, ...).
Running the scripts rebuilds every figure and table.

## R environment

R 4.4.2 on Ronen's machine. Packages are not pinned with `renv`.

CRAN: `tidyverse`, `magrittr`, `vegan`, `ape`, `phangorn`, `writexl`, `readxl`,
`cowplot`, `patchwork`, `ggh4x`, `flextable`, `officer`, `dunn.test`,
`ggsignif`, `RColorBrewer`, `magick`, `httr`, `xml2`, `jsonlite`, `here`,
`microeco`, `file2meco`.
Bioconductor: `Biostrings`, `ggtree`, `ANCOMBC`, `phyloseq`.
GitHub: `qiime2R` (`jbisanz/qiime2R`).

The BLAST scripts send an e-mail address to NCBI, as NCBI asks. Set it before
running them: `Sys.setenv(NCBI_EMAIL = "you@example.org")`. The PowerShell
runners look for MAFFT, IQ-TREE and BLAST+ in `$HOME\tools` (override with
`$env:OCTO_TOOLS`).

## Conventions in the code

`set.seed(123)` before any random process. `clone()` before subsetting a
microtable object. `%<>%` for in-place modification. `# ===` section comments.
The species palette is defined once, in `00_setup.R` (`species_colors`), and
matches Figure 2.
