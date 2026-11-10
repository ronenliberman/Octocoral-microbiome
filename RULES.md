# Processing rules

A living record of the data-handling and analytical decisions in this project,
with the reason for each, so that a decision can be revisited on the record
rather than silently. Update this file in the same commit as the code change it
describes. The full reasoning and the review history are in
`documentation/PROJECT_CONTEXT.md`.

Every rule below is what the published paper (FEMS Microbiology Ecology,
fiag124) used. Changing one changes a published result (AGENTS.md §3).

---

## Global rules

### Input data
- **One ASV table for everything: PEE2.** `raw files/table-NoEMC-PEE2.qza`
  (DADA2 feature table, chloroplast and mitochondria removed),
  `taxonomy-PEE2.qza` (SILVA naive Bayes classifier, 515F/806R),
  `dna-sequences-PEE2.fasta` (representative sequences), and
  `metadata_all_edited.csv`. Imported with `file2meco::qiime2meco()`.
  The R1 trees and the R1 *Endozoicomonas* master table were built from an older
  table (`table-NoEMC.qza`); round 3 moved everything to PEE2, and the older
  table is no longer in the repository (it is in the git history).
- `Figure4_species_ANCOMBC_REVISION2.R` still names `raw data/table-NoEMC.qza`
  and two exclusions. It is never run directly: `06c` substitutes the PEE2
  files and the three exclusions. Re-running ANCOM-BC2 on PEE2 (1 Oct 2026) gave
  the same 24 year and 9 site results, with identical log fold changes.

### Sample filters (same in every script)
- Kingdom restricted to Bacteria and Archaea.
- Hosts `BAST`, `ECAR`, `EFLE`, `MMUR` (plus `SW` where seawater is analysed).
- Site `DC2` is recoded to `DC3`; sites kept: `DC1`, `DC3`, `BC4`, `BC6`.
- **Three samples are excluded:** `ECAR_DC3_3_2_2023`, `BAST_DC3_1_2024`,
  `EFLE_DC3_3_2024` (below 20,000 reads).
- Result: 180 octocoral samples, 42,576 ASVs.

### Normalization, by test
| Test | Data |
|------|------|
| Alpha diversity (richness, evenness) | Rarefied to **20,000 reads/sample** |
| Community composition, core microbiome | TSS relative abundance |
| Beta diversity, PCoA, PERMANOVA | robust CLR (`rclr`), Aitchison distance |
| ANCOM-BC2 | Unrarefied filtered counts |

### Statistics
- PERMANOVA with `vegan::adonis2`, sequential terms. The script uses 999
  permutations; see PROJECT_CONTEXT.md "Open items" for the Table 2 footnote.
- ANCOM-BC2 within each host for year (reference 2023) and for site; FDR,
  alpha = 0.05.
- `set.seed(123)` before any random process.

### Species palette
`B. asbestinum` `#E41A1C`, `E. caribaeorum` `#FF7F00`, `E. flexuosa` `#984EA3`,
`M. muricata` `#4DAF4A`, seawater `#377EB8`. Defined once in
`script/Editor_round3/00_setup.R`. The R1 Figure 3 had two colours swapped;
every round-3 figure uses this vector.

### ASV labels
- **One stable label per ASV across every figure and table** (`02` writes
  `ASV_stable_IDs.csv`; it stops if a label is duplicated).
- The six core *Endozoicomonas* ASVs keep their R1 labels Endo_1–Endo_6, fixed
  by sequence hash in `00_setup.R` (`core_endo_map`).
- Myco_1–Myco_10 are the ten most abundant *Mycoplasma* ASVs by total reads,
  taken from the R1 master table (`blast/Mycoplasma_master_table.csv`).
  Myco_11 and higher were renumbered on PEE2, so R1 numbers above 10 do not
  carry over.
- **Style:** the main text writes "Endo 1"; figures, the supplement and the
  workbook write "Endo_1". This is deliberate.

---

## ASV verification (Endozoicomonadaceae, Mycoplasmataceae)

**Scripts:** `01`, `01b`, `01b_Myco`, `01c`, `01d`, `02`.

Thresholds from Yarza et al. (2014, *Nat Rev Microbiol* 12:635): median 16S
identity within a family **86.5%** and within a genus **94.5%**. They come from
near-full-length 16S, so for 253 bp V4 reads they are approximate.

| Status | Rule (best hit in the NCBI 16S type-strain database) |
|---|---|
| Verified (genus) | SILVA family = NCBI family, identity ≥ 94.5%, SILVA genus = NCBI genus |
| Verified (family) | SILVA family = NCBI family, identity ≥ 86.5% |
| Unverified | anything else; excluded from the trees and from genus counts |

- **Tier B.** Unverified ASVs are searched against `core_nt` (which has the
  uncultured coral clones the type-strain database lacks) and rescued only if
  their best hit that has a family-level lineage passes the same rule.
- **Genus check.** Family-level ASVs are searched against `core_nt` again and
  upgraded to genus only on the same rule. Twelve family-level ASVs were never
  tested at genus level because the NCBI website returned partial results; they
  keep the family label.
- **Family equivalents.** SILVA 138 uses the pre-2018 Mollicutes taxonomy, so
  `Metamycoplasmataceae` and `Mycoplasmoidaceae` count as `Mycoplasmataceae`.
- **Class labels for unverified ASVs.** An unverified ASV whose best `core_nt`
  hit is ≥ 78.5% (the class median) is labelled by that class, for example
  "c_Spirochaetia (Unv_1)".
- **Mycoplasma host assignment** counts a colony only when *Mycoplasma* is more
  than 1% of its reads (`MYCO_MIN_FRACTION`).
- The searches are saved in `blast/` and reused. Re-querying NCBI would give
  different hits as its databases grow.

---

## Bacterial trees (Figs 5–6)

**Scripts:** `03`, `04`, `11`, `16`, `16b`; helpers in `00_setup.R`.

- **Tips:** verified ASVs detected in ≥ 5 samples, plus all core ASVs, plus
  references.
- **References are RefSeq type-strain records (NR_) only,** looked up by
  organism name and listed in `00_setup.R` (`endo_references`,
  `myco_references`). The script refuses a record whose NCBI title does not
  match its label. Seven of the fourteen R1 reference accessions were not the
  labelled species, which is why they were replaced.
- References are trimmed to the 515F/806R amplicon so they match the ASVs.
- Outgroups: *Kistimonas asteriae* (Endozoicomonadaceae),
  *Acholeplasma laidlawii* (Mycoplasmataceae).
- MAFFT L-INS-i (`--localpair --maxiterate 1000`); IQ-TREE ModelFinder
  (`-m MFP`), 1,000 ultrafast bootstraps and 1,000 SH-aLRT replicates.
- **Support shown as SH-aLRT/UFBoot** at nodes with UFBoot ≥ 70%; a node counts
  as supported at SH-aLRT ≥ 80% and UFBoot ≥ 95% (Guindon et al. 2010; Minh et
  al. 2013). The trees are rooted with `edgelabel = TRUE` so labels stay on
  their nodes.

---

## Core microbiome (Fig 3, Appendix 2)

**Scripts:** `Core_Microbiome_Revised_AllThresholds.R`, `05`, `08`.

- Prevalence thresholds 80, 90 and 100% of colonies per host, on TSS data.
- Figure 3 shows only core ASVs with ≥ **1%** mean relative abundance within
  the host (`MIN_MEAN_RELABUND_FIG3`); Appendix 2 lists all of them.
- The 90% core is checked by 1,000 random subsamples of 20 colonies per host,
  **without replacement**, seed 123 (`08`). An ASV is robust if it is core in
  ≥ 80% of subsamples; 70–80% is reported as borderline.

---

## Reopened or pending decisions

None reopened. See PROJECT_CONTEXT.md "Open items".
