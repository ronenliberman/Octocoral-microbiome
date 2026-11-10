# Octocoral Microbiome — Manuscript Revision Analysis Knowledge Base

## Purpose

This document captures the complete analytical pipeline, decisions, and reasoning behind the revision of:

> **"Host species identity drives temporally stable octocoral microbiomes in Florida"**
> Manuscript ID: **FEMSEC-25-12-0321** (FEMS Microbiology Ecology)
> Authors: Ronen Liberman (RL), Matías Gómez-Corrales (MGC), Jose V. Lopez (JVL)

It is intended as a standing reference for the associated GitHub repository and for anyone working on this codebase later, so that analytical choices made during peer review do not have to be reconstructed from scratch. It covers the first revision (FEMSEC-25-12-0321) and predates the round-3 changes recorded in `documentation/PROJECT_CONTEXT.md`.

---

## 1. Study System & Data

- **Host species (4 Florida octocorals):**
  - *Briareum asbestinum* (BAST) — encrusting
  - *Erythropodium caribaeorum* (ECAR) — encrusting
  - *Eunicea flexuosa* (EFLE) — branching
  - *Muricea muricata* (MMUR) — branching
- **Sites:** BC4, BC6, DC1, DC3 (Broward County / Deerfield Beach, Florida)
- **Years:** 2023, 2024, 2025 (plus seawater controls, SW)
- **Sequencing:** Illumina MiSeq 2×300 bp, 16S rRNA V4 region, EMP 515F/806R primers, SILVA 138.2 reference
- **Consistent species color scheme (used across all figures):**

```r
species_colors <- c(
  "B. asbestinum"  = "#E41A1C",  # BAST
  "E. caribaeorum" = "#FF7F00",  # ECAR
  "E. flexuosa"    = "#984EA3",  # EFLE
  "M. muricata"    = "#4DAF4A"   # MMUR
)
```

---

## 2. Reviewer Comments — What Drove the Revision

Two reviewers returned major comments. The revision was built directly around resolving these:

| Issue | Reviewer ask | Resolution |
|---|---|---|
| Dominant phyla reported without stats | Add statistical comparisons; move to family-level | Family-level ANCOM-BC2 + explicit stats added |
| Core microbiome results selective | Only 90% threshold shown despite testing multiple | All three thresholds (80/90/100%) now reported, with bootstrap validation |
| ANCOM-BC2 vague ("several," "multiple") | Report explicit counts, effect sizes, ASV-level results | Re-ran ANCOM-BC2 at ASV level on unrarefied data; explicit LFC/q-value tables |
| Morphology effects discussed but not tested | Show results, not just a supplementary figure | PERMANOVA confirmed morphology is *nested within* species (not an independent factor — see Section 5) |
| Inconsistent taxonomic levels (ASV vs. family vs. genus) | Be consistent | Core microbiome and ANCOM-BC2 both fixed at ASV level throughout |
| *Mycoplasma* discussed but no results shown | Add explicit results | Added dedicated *Mycoplasma* phylogenetic/network analysis |
| Justify why only 90% shown | Justify or show all | All thresholds shown; 90% presented as primary because of bootstrap validation |

---

## 3. Bioinformatics Pipeline (QIIME2 → R)

### 3.1 QIIME2 processing
- Platform: QIIME2 2023.2 (or matching version across merged runs)
- Reads imported via manifest, `dada2 denoise-paired`
- **Revised DADA2 parameters** (vs. original submission):
  - `--p-max-ee 2` (tightened from an original essentially-unfiltered `--p-max-ee 29`)
  - `--p-trunc-len-f` / `--p-trunc-len-r` explicitly specified per run (quality-plot-guided)
  - Chimera removal via consensus method
- Taxonomy: `classify-sklearn` with SILVA 138.2 (515F/806R-trained classifier)
- Filtering: exclude Eukaryota, Mitochondria, Chloroplast (`qiime taxa filter-table`)
- Multi-year runs merged with `qiime feature-table merge` / `merge-seqs` (only valid when using identical QIIME2 version and trimming parameters across runs)
- Exported to `.tsv`/`.fasta` via `qiime tools export` + `biom convert`

### 3.2 Import into R (microeco via qiime2R/file2meco)

```r
data <- qiime2meco(
  "raw_data/table-NoEMC-PEE2.qza",
  taxonomy_table = "raw_data/taxonomy-PEE2.qza",
  sample_table   = "raw_data/metadata_all_edited.csv",
  auto_tidy      = TRUE
)

# Restrict to prokaryotes
data$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
```

### 3.3 Sample/site filtering and cleanup

```r
d <- clone(data)
d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE", "SW"))
d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)   # site code correction
d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

# Three samples removed for the revision due to unreliable sequencing depth/quality
samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024", "EFLE_DC3_3_2024")
d$sample_table <- d$sample_table %>% filter(!rownames(.) %in% samples_to_remove)

species.names <- list(
  "SW"   = "Seawater",
  "BAST" = "B. asbestinum",
  "ECAR" = "E. caribaeorum",
  "MMUR" = "M. muricata",
  "EFLE" = "E. flexuosa"
)
d$sample_table <- d$sample_table %>% mutate(SPP = recode(SPP, !!!species.names))
d$tidy_dataset()
```

**Note:** these same filtering steps (species/site subset, Site DC2→DC3 correction, sample removal, species renaming) are repeated identically at the top of every downstream script (core microbiome, ANCOM-BC2, ordinations, phylogenetics) to guarantee a consistent dataset across all analyses.

---

## 4. Normalization Strategy (Revision Change — Reviewer-Driven)

The single biggest methodological shift in the revision was moving from a single Bray-Curtis-on-rarefied-data approach to **method-matched normalization for each analysis type**:

| Purpose | Method | Rationale |
|---|---|---|
| Alpha diversity | Rarefaction to 20,000 reads/sample | Standard for richness/diversity estimators; depth must be equalized |
| Community composition / relative abundance plots | TSS (total sum scaling) | Avoids discarding reads via rarefaction when only proportions are needed |
| Beta diversity | rCLR (robust centered log-ratio) → Aitchison distance (Euclidean on rCLR) | Compositionally appropriate; avoids Bray-Curtis's known distortions with compositional data; replaces original Bray-Curtis approach |
| Differential abundance (ANCOM-BC2) | Raw, unrarefied counts | ANCOM-BC2 performs its own bias correction internally; rarefying first would discard information and is not recommended for this method |

```r
mt_rarefied <- trans_norm$new(dataset = d)$norm(method = "rarefy", sample.size = 20000)
mt_TSS      <- trans_norm$new(dataset = d)$norm(method = "TSS")
mt_CLR      <- trans_norm$new(dataset = d)$norm(method = "rclr")
```

For beta diversity, Euclidean distance is computed on the rCLR-transformed table (this *is* the Aitchison distance):

```r
euc_dist <- vegan::vegdist(t(mt_CLR$otu_table), method = "euclidean")
mt_CLR$beta_diversity <- list(euclidean = as.matrix(euc_dist))

t_beta_clr <- trans_beta$new(dataset = mt_CLR, group = "SPP", measure = "euclidean")
t_beta_clr$cal_ordination(method = "PCoA")
```

Rarefied Bray-Curtis PERMANOVA was retained *alongside* the Aitchison-based analysis (not replaced) so both are reported for transparency/comparison with prior microbiome literature.

---

## 5. PERMANOVA / Beta Diversity Statistics

- Hierarchical design: overall species effect first, then nested year/site effects **within** each species.
- PERMDISP (999 permutations) run first to confirm homogeneity of dispersion before interpreting PERMANOVA.
- Pairwise comparisons (species × species, species × seawater) use Bonferroni correction.
- Marginal sums of squares used for multi-factor models (`manova_all = TRUE` equivalent / `by = "margin"` in `adonis2`).

**Key learning (important for interpretation):** Morphology is **perfectly nested within species** for this dataset (each species has one morphology type: BAST/ECAR = encrusting, EFLE/MMUR = branching). When morphology is included as a factor alongside species, it explains **R² = 0.000 (marginal)** — i.e., it adds no independent explanatory power once species is in the model. This is *not* a bug; it is the correct statistical outcome and directly answers the reviewer's morphology question: morphology effects observed in earlier exploratory analyses were actually species effects. Do not treat morphology as an independent covariate in models that already include species.

Effect size interpretation used throughout: R² > 0.25 large, 0.09–0.25 medium, < 0.09 small (Kelly et al. 2015 convention).

---

## 6. Core Microbiome Analysis

### 6.1 Multi-threshold prevalence approach
Following Van Oppen & Blackall (2019), core taxa defined at the ASV level using **three prevalence thresholds per species**, reported side by side (this directly answers the reviewer complaint that only 90% was originally shown):

- **80%** (moderate)
- **90%** (stringent — primary reported threshold)
- **100%** (sensu stricto)

Prevalence = % of samples within a species containing the ASV (presence/absence). Mean relative abundance = average proportional abundance of the ASV across samples within that species.

### 6.2 Bootstrap validation of the 90% core (key methodological addition)

Rationale: a raw 90% prevalence threshold on the full sample set can be sensitive to sampling variability. Bootstrap resampling tests whether an ASV remains "core" under repeated random subsampling.

```r
n_bootstrap <- 1000
n_samples_per_bootstrap <- 20     # sampled without replacement, per species, per iteration
bootstrap_threshold_cutoff <- 80  # % of iterations an ASV must appear as core to be called "robust"

bootstrap_core_microbiome <- function(microtable_obj, prevalence_threshold,
                                       n_iterations, n_samples_bootstrap) {
  all_samples <- colnames(microtable_obj$otu_table)
  bootstrap_core_lists <- list()

  for (i in 1:n_iterations) {
    sampled_samples <- sample(all_samples, n_samples_bootstrap, replace = FALSE)
    temp <- clone(microtable_obj)
    temp$sample_table <- temp$sample_table[sampled_samples, , drop = FALSE]
    temp$otu_table     <- temp$otu_table[, sampled_samples, drop = FALSE]
    temp$tidy_dataset()

    core_result <- calculate_prevalence_core(temp, prevalence_threshold)
    bootstrap_core_lists[[i]] <- core_result$core_taxa
  }

  all_asvs <- unique(unlist(bootstrap_core_lists))
  asv_frequency <- sapply(all_asvs, function(asv) {
    sum(sapply(bootstrap_core_lists, function(x) asv %in% x))
  })
  data.frame(
    ASV_ID = names(asv_frequency),
    Bootstrap_frequency_percent = (asv_frequency / n_iterations) * 100
  )
}
```

An ASV is called part of the **"robust core"** only if it meets the 90% prevalence threshold *and* appears as core in ≥80% of the 1,000 bootstrap iterations. Bootstrap frequency is reported alongside prevalence and mean relative abundance for every robust core ASV (see `Robust_Core_90pct_Bootstrap_Summary.csv`).

### 6.3 Output structure
- `Core_Microbiome_All_Thresholds_All_Species.xlsx` — one sheet per species × threshold combination
- `Robust_Core_90pct_Bootstrap_Summary.csv` — final validated core, columns: Species, ASV_ID, Genus_ID (globally numbered per genus, e.g. "Endozoicomonas 1"), Class, Family, Genus, Prevalence, Mean_Relative_Abundance_Percent, Relative_Abundance_Total_Microbiome, Bootstrap_frequency_percent

### 6.4 Global Genus ID numbering (for readability in figures/tables)
When multiple ASVs from the same genus appear in the robust core, they are numbered sequentially (Genus 1, Genus 2, …) so figures can refer to "Endozoicomonas 1–6" rather than raw ASV hashes:

```r
genus_numbering <- all_robust_core %>%
  select(ASV_ID, Genus) %>% distinct() %>%
  arrange(Genus, ASV_ID) %>%
  group_by(Genus) %>%
  mutate(Genus_Number = row_number(), Genus_ID = paste0(Genus, " ", Genus_Number)) %>%
  ungroup() %>% select(ASV_ID, Genus_ID)
```

---

## 7. Differential Abundance — ANCOM-BC2

- Package: `ANCOMBC` (Lin & Peddada 2020)
- Run on **raw, unrarefied ASV-level counts** (bias correction is internal to the method — do not pre-rarefy)
- FDR (BH) correction across ASVs; significance defined as q < 0.05
- **Output structure is long-format, not wide:** single `lfc`, `p_val`/`P.adj`, `diff` columns with a `Factors` column indicating which pairwise comparison each row belongs to. Downstream code (`res_diff` object) must be filtered/reshaped accordingly rather than assuming one column per comparison.
- Explicit counts of significant ASVs per comparison, plus log-fold changes, are now reported in a supplementary table (directly resolving the reviewer's request for quantitative rigor instead of "several"/"multiple").

---

## 8. Phylogenetic Analyses of Dominant Symbionts

Two genera were selected for targeted phylogenetic + haplotype network analysis because of their dominance in the core microbiome: ***Endozoicomonas*** and ***Mycoplasma***.

### 8.1 Common pipeline
- Extract all ASVs assigned to the genus (and family, to catch unclassified-at-genus ASVs) from the filtered dataset
- Align representative sequences (MAFFT L-INS-i)
- V4 region trimming via regex anchors: forward `CAGC.GCCGCGGTAA`, reverse `GTAGTCC`
- Tree building: Neighbor-Joining + Jukes-Cantor (JC69) model + 1000 bootstrap replicates (`phangorn`)
- Haplotype networks: `pegas` + `ggplot2`/`ggforce`
- NCBI BLAST validation via web interface, Hit Table CSV export (no header row)
- BLAST hits ranked by bitscore; top 3 hits retained per ASV with an `identity_category` annotation rather than a hard 97% cutoff (many hits fall below 97% and are still informative)

### 8.2 *Endozoicomonas*
- 90 ASVs extracted from *Endozoicomonas*/Endozoicomonadaceae
- 6 core ASVs labeled Endo_1–Endo_6 for the manuscript figures
- BLAST identifications:
  - Endo_2 = *E. euniceicola* (100%)
  - Endo_4, Endo_5 = *E. euniceicola* (99.6%)
  - Endo_3, Endo_6 = *E. lisbonensis* (98.8–99.6%)
  - Endo_1 = *E. atrinae* (98.8%)
- Two clear lineages recovered: *E. euniceicola*-like and *E. gorgoniicola*-like
- 24/90 ASVs shared across more than one host species; 38/90 shared across more than one site — indicating this genus is not strictly host-specific at the strain level despite being core in every species
- All 10 currently described *Endozoicomonas* species represented among BLAST hits
- **No 1% abundance filtering applied** — dominant across all four hosts, so trace-occurrence noise is not a major concern
- Outgroup: *Kistimonas*

### 8.3 *Mycoplasma*
- Dominant in the two branching species: EFLE (~34% relative abundance), MMUR (~29%)
- Only trace abundance in the encrusting species: BAST (~0.19%), ECAR (~0.07%)
- All recovered ASVs are divergent from described type strains (*M. iowae*, *M. insons*, *M. todarodis*) — consistent with novel, octocoral-specific lineages rather than known pathogenic/commensal *Mycoplasma* species
- Two clades recovered: an upper clade dominated by MMUR, and a lower, mixed MMUR+EFLE clade
- **1% relative abundance threshold applied before tree tip coloring and haplotype network construction** — this is a deliberate, non-default choice: including trace occurrences from the encrusting species (which are almost certainly incidental/background reads rather than colonizing symbionts) would add biologically uninformative noise to the tree and network. This is the opposite decision from the *Endozoicomonas* analysis and should not be applied mechanically across genera — the filtering choice depends on whether trace abundance in non-dominant hosts is biologically meaningful.
- Outgroup: *Acholeplasma*

---

## 9. Working Environment & File Organization

- **Two-machine workflow:**
  - Mac terminal: QIIME2, MAFFT, IQ-TREE, alignment steps
  - Windows PC / RStudio: all R analysis and visualization
  - Files transferred manually between machines as needed
- **Working directory (revision):** `~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/`
  - Raw data: `raw_data/` (QIIME2 `.qza` exports, `metadata_all_edited.csv`)
  - Outputs: `Output/` subdirectory
- **Key final scripts:**
  - `Octocoral_Revised_Analysis_COMPLETE.R` — full revised pipeline (import → normalization → alpha/beta diversity → PERMANOVA)
  - `Core_Microbiome_Revised_AllThresholds.R` — multi-threshold core microbiome (no bootstrap)
  - `Core_Microbiome__90__Prevalence_with_BOOTSTRAPING.R` — bootstrap-validated 90% core
  - `Figure4_species_ANCOMBC.R` — ANCOM-BC2 analysis and figure
  - `species-specific_PCoA_ordinations_Nov2025.R` — per-species rarefied PCoA with sample labels, colored by Year, shaped by Site
  - `Endozoicomonas_Phylogenetic_Analysis_FINAL.R` (also `_CLEAN.R`) — Endozoicomonas tree + haplotype network
  - `Mycoplasma_Phylogenetic_Analysis.R` — Mycoplasma tree + haplotype network
- **Documents produced for submission:**
  - `Octocoral_Track_Changes.docx` (384 paragraphs) — Word track-changes version
  - Rebuttal letter — comment-and-response format (black = reviewer text, red = author responses)
  - Cover letter
  - Preprint deposited on bioRxiv (CC BY 4.0)
  - Raw sequence data: NCBI BioProject **PRJNA1373357**

---

## 10. R Coding Conventions Used Throughout

- `set.seed(123)` before any randomization step (bootstrapping, ordination with random starts, etc.)
- `dplyr::select()` / `dplyr::first()` explicitly namespaced to avoid masking conflicts (esp. with `MASS`, `Biostrings`, `phyloseq`)
- Sequential `cat()` statements as lightweight progress logging through long scripts, organized under `SECTION`/`STEP`/`PHASE` headers
- Bioconductor S4 objects (e.g., from `phyloseq`, `Biostrings`) converted with `as.data.frame()` and stripped of list-columns before piping into `dplyr`
- `otu_table` single-row subsetting returns a `data.frame` rather than a named vector in this codebase's usage pattern — always wrap with `unlist()` when a numeric vector is needed
- IQ-TREE modifies tip labels (spaces become underscores; long multi-host names get truncated) — tip-to-metadata joins use an explicit `for` loop rather than `sapply`/`merge` to avoid silent mismatches; node support labels in `SH-aLRT/UFboot` format are split on `/` to extract UFboot separately
- QIIME2 metadata/ID text files require the header `feature-id` for ASV ID columns
- `.docx` files in this project's working directories are sometimes plain-text markdown saved with a `.docx` extension (not true OOXML) — `python-docx`/`zipfile`-based readers will fail on these; convert with `pandoc` first, then parse as markdown

---

## 11. Writing & Figure Conventions

- Results sections use results-oriented language only — no embedded methods or discussion
- Every interpretive claim carries a citation
- Temperate language throughout ("suggest," "indicate") rather than "demonstrate"/"prove"
- No em-dashes; American English spelling
- Figures: 300 dpi, consistent species color scheme (Section 1), legends positioned to avoid obscuring data points, x-axis labels rotated 45° where category labels risk overlapping

---

## 12. Key References (Revision-Relevant)

- Callahan et al. (2016) — DADA2
- Lin & Peddada (2020) — ANCOM-BC2
- Liu et al. (2021) — microeco R package
- Bolyen et al. (2019) — QIIME2
- Van Oppen & Blackall (2019) — core microbiome definitional framework
- Kelly et al. (2015) — effect size interpretation convention for R² in community ecology

---

## 13. Status Snapshot (as of this revision)

- Revision substantially complete: all reviewer major comments addressed with corresponding methodological or reporting changes (Section 2 table)
- Remaining step: final submission of revised manuscript, rebuttal letter, and supplementary materials to FEMS Microbiology Ecology
- This dataset (3 years, 4 species) also now serves as the **baseline/comparison dataset** for a second, analytically distinct project on bleaching-associated microbiome shifts, and will be referenced (not re-analyzed) in an upcoming ICRS 2026 conference talk
