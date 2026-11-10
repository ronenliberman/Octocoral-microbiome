# Project context

The science record for this repository: what the study is, how the analysis
got to its published form, which decisions are settled, and what is still open.
`RULES.md` is the short index of the decisions; this file says why.

## 1. The paper

Liberman R, Gómez-Corrales M, Lopez JV. Host species identity drives temporally
stable octocoral microbiomes in Florida. *FEMS Microbiology Ecology*, fiag124
(manuscript FEMSEC-26-04-0109; an earlier submission was FEMSEC-25-12-0321).
Accepted in 2026; proofs corrected in November 2026.

## 2. Design

- **Hosts:** *Briareum asbestinum* (BAST) and *Erythropodium caribaeorum* (ECAR),
  encrusting; *Eunicea flexuosa* (EFLE) and *Muricea muricata* (MMUR),
  branching.
- **Sites:** four reef sites off southeast Florida, a shallow and a deep site in
  Broward and in Miami (codes `BC4`, `BC6`, `DC1`, `DC3`; coordinates and depths
  in Table S1).
- **Years:** 2023, 2024, 2025; seawater controls (`SW`).
- **Sequencing:** 16S rRNA V4 (515F/806R), Illumina MiSeq; DADA2 in QIIME2;
  SILVA taxonomy. Reads in BioProject PRJNA1373357.
- **Samples analysed:** 180 octocoral colonies after the three exclusions in
  `RULES.md`.

## 3. Review history

**Round 1** (FEMSEC-25-12-0321) asked for consistent normalization, ASV-level
reporting, NCBI validation of poorly classified ASVs, Aitchison distances and
host-specific differential abundance. These shaped the round-2 scripts in
`script/` (March–August 2026). `Octocoral_Manuscript_Revision_Analysis_Knowledge.md`
in this folder is the record of that round.

**Round 2** (FEMSEC-26-04-0109) added colour-blind-safe Figure 4 styling, a
bubble version of Figure 4 (the one published) and Table S8 on base-pair
differences among core ASVs.

**Round 3** (editor's minor revision, September–October 2026) is the reason for
`script/Editor_round3/`. The editor's points and what they changed:

| Point | What was wrong | What changed |
|---|---|---|
| Tree clades unlabelled | Only the core ASVs had labels | Every tip labelled with a stable ID |
| Tree squeezed to one side | Ten "*Endozoicomonas*" ASVs were only 75–81% identical to any type strain, and 7 of the 14 R1 reference accessions were not the labelled species (for example "*E. lisbonensis*" was a *Vibrio*, the *Mycoplasma* outgroup was *Bordetella*) | Unverified ASVs removed by the Yarza-threshold rule; references replaced by RefSeq type strains checked by name; trees rebuilt |
| Haplotype network repeats the tree | | Network panels and text removed |
| Classification at 75% identity suspect | Also an R1 rank bug in Appendix 1 (the family column held the order) | Two-tier verification (type strains, then `core_nt`); Appendix 1 rebuilt with ranks looked up by name |
| The same ASV numbered differently in Figs 3, 4 and 5 | Each figure numbered its own ASVs | One stable label per ASV (`02`) |
| Aligned positions not stated | | Alignment length, gap-free positions and parsimony-informative sites reported |
| Figure 3 mostly empty | All ≥ 80% prevalence ASVs drawn | Only core ASVs with ≥ 1% mean abundance |
| Supplementary tables not received | Upload designation | One workbook (`08`) |

While re-examining the data in round 3:

- ANCOM-BC2 was re-run on the PEE2 table with the three standard exclusions;
  the results were identical to round 2 (24 year and 9 site results).
- The 90% core bootstrap gave 22/29 robust ASVs in *E. caribaeorum* (the
  earlier text said 24); seven ASVs at 91.1% prevalence score 75.6–78%, just
  under the 80% cut-off, and are reported as borderline.
- Tree support was first shown as UFBoot only with classical bootstrap cut-offs,
  which overstates it; the published trees show SH-aLRT/UFBoot with the
  SH-aLRT ≥ 80% and UFBoot ≥ 95% criterion. A rooting bug (`edgelabel`) that
  misplaced support labels was fixed before submission.

## 4. Results the scripts reproduce

| Result | Script | Output |
|---|---|---|
| Alignments: Endozoicomonadaceae 33 sequences (23 ASVs + 9 type strains + *Kistimonas asteriae*), 255 sites, 251 gap-free, 28 parsimony-informative, HKY+R2; Mycoplasmataceae 15 sequences (11 ASVs + 3 type strains + *Acholeplasma laidlawii*), 255 sites, 243 gap-free, 78 parsimony-informative, TVM+I | `03`, `04` | `Alignment_stats_*.csv` |
| Pairwise identity among core ASVs: Endo_1–6 95.6–99.6%; Myco_1–10 81.2–99.6% | `07` | Table S8 sheet of the workbook |
| Verification of all 131 Endozoicomonadaceae and Mycoplasmataceae ASVs | `01`, `07` | Table S9 sheet, Appendix 1 |
| Robust 90% core: *B. asbestinum* 6/7, *E. caribaeorum* 22/29, *E. flexuosa* 1/1, *M. muricata* 2/2 | `08` | Table S10 sheet |
| Figure 3: 17 rows (6, 7, 3 and 1 per host) | `05` | Figure 3, Appendix 2 |

The workbook sheet names follow the round-3 numbering (S8–S10); the published
supplement renumbered the tables.

## 5. Settled decisions

See `RULES.md`. The ones most likely to be questioned again:

- The PEE2 table for every analysis.
- The Yarza et al. (2014) thresholds, applied to V4 reads as approximations.
- Endo_1–6 and Myco_1–10 keep their R1 labels; other numbers were reassigned.
- "Endo 1" in the main text, "Endo_1" everywhere else.
- Endo_6 sits next to *E. gorgoniicola* in the tree on an unsupported node; by
  sequence it is closest to *E. lisbonensis*. The paper says only that Endo 3
  and Endo 6 fall in a poorly resolved part of the tree.

## 6. Open items

- **Table 2 permutations.** `Octocoral_Revised_Analysis_COMPLETE.R` runs
  PERMANOVA with 999 permutations. The published Table 2 footnote says 9,999.
  With 9,999 the p-values would be: Species 0.0001; Year 0.0002–0.0004 (\*\*\*
  instead of \*\*); Site 0.003–0.005; Species × Site 0.0001–0.0002 (\*\*\*
  instead of \*\*). Decide whether the script or the footnote is the record.
- `17_Fig2_y_axis_title.R` starts from a GIMP-assembled Figure 2, so Figure 2
  cannot be rebuilt from code alone.
- `Supplementary_Materials.R` and the community-composition script need objects
  from `Octocoral_Revised_Analysis_COMPLETE.R` in memory; they are not
  standalone.
- Package versions are not pinned (`renv`).
