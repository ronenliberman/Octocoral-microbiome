# phylogeny/bacteria

The two bacterial trees of the paper: Endozoicomonadaceae (Figure 5, prefix
`Endo_r3`) and Mycoplasmataceae (Figure 6, prefix `Myco_r3`). Settings and
reference selection are in RULES.md, "Bacterial trees".

| File | What it is | Written by |
|---|---|---|
| `*_for_tree.fasta` | Tree input: verified ASVs (≥ 5 samples, plus core ASVs) and RefSeq type-strain references, trimmed to the V4 amplicon | `03` / `04`, phase 1 |
| `*_reference_names.csv` | Reference accession, label and outgroup flag | `03` / `04`, phase 1 |
| `*_mafft_aligned.fasta` | MAFFT L-INS-i alignment | MAFFT |
| `*_iqtree.treefile` | Maximum-likelihood tree; node labels are SH-aLRT/UFBoot | IQ-TREE |
| `*_iqtree.contree` | UFBoot consensus tree | IQ-TREE |
| `*_iqtree.iqtree`, `*_iqtree.log` | IQ-TREE report (ModelFinder, best model) and log | IQ-TREE |

`11_run_alignment_and_tree.ps1` reruns MAFFT and IQ-TREE on Windows. The trees
here were inferred on 1 Oct 2026 with IQ-TREE 3.1.4 (`--seed 123`); the
command is at the top of each `.log`.

Re-running `03`/`04` rewrites the two phase-1 files; `git diff` should show no
change. Re-running MAFFT or IQ-TREE replaces the published trees, so do that
only through an Analysis issue (AGENTS.md §3). Figures 5 and 6 are drawn from
these files by `16_tree_figures_v3.R` and `16b_Fig5_v3_legend_positions.R`.
