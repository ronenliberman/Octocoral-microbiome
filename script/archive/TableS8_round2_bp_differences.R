# Archived from script/Supplementary_Materials.R (round 2, August 2026).
# Superseded by script/Editor_round3/07_TableS8_revised.R, which computes the
# published Table S8 from the round-3 alignments. Paths are the original ones.

################################################################################
# Supplementary Table S8: Pairwise base pair differences AMONG core
# Endozoicomonas (Endo_1-6) and Mycoplasma (Myco_1-10) ASVs.
#
# Author: Ronen Liberman
# Date: August 2026
# Purpose: Addresses Editor comment (FEMSEC-26-04-0109, Round 2):
#   "it would help to state the number of base pair differences between
#   Endo 2-6 and other taxa in the supplemental material"
#   (clarified: differences are wanted among the core ASVs themselves,
#   e.g. Endo_1 vs. Endo_2, Endo_1 vs. Endo_3, etc. - not vs. reference
#   type strains)
#
# Dependencies: This script REUSES objects already built by
#   Endozoicomonas_Phylogenetic_Analysis_FINAL.R and
#   Mycoplasma_Phylogenetic_Analysis.R. It does not repeat BLAST, ASV
#   extraction, or the master-table construction - it reads the same
#   aligned FASTA files those scripts used to build Figures 5A/6A, and
#   only performs the new pairwise-distance calculation and table
#   export requested here.
#
# Output: Supplementary_Table_S8_bp_differences.docx, containing two
#   tables:
#     1. Pairwise bp differences among the 6 core Endozoicomonas ASVs,
#        and among the 10 core Mycoplasma ASVs (full pairwise matrix,
#        long format)
#     2. Same pairwise comparisons summarized against the full ASV set
#        in each tree (i.e., how divergent is each core ASV from the
#        full range of other ASVs in its genus, not just other core ASVs)
################################################################################

# ---- PACKAGES ----------------------------------------------------------------
library(ape)
library(Biostrings)
library(tidyverse)
library(officer)
library(flextable)

# ---- PATHS (reuse existing project conventions) -------------------------------
# These match the paths already defined in the Endo/Myco phylogenetic scripts.
base_dir    <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/16 rRNA analysis/"
rev_dir     <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/"
output_dir  <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Output/"
phylo_dir   <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/raw_data/phylo analysis"
endo_output <- file.path(output_dir, "Endozoicomonas_phylogeny")
myco_output <- file.path(output_dir, "Mycoplasma_phylogeny")

supp_output <- file.path(output_dir, "Supplementary_Tables")
dir.create(supp_output, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1: Load the SAME aligned sequences used to build Figures 5A and 6A
# (read from phylo_dir, matching what tree_bs / tree_bs_myco were built
# from in the Endo/Myco scripts — do NOT re-align, to keep bp differences
# consistent with what the published trees actually show)
# =============================================================================

cat("\n========== STEP 1: Loading tree alignments ==========\n\n")

# --- Endozoicomonas: MAFFT-aligned FASTA (aligned on Mac, per project workflow) ---
endo_aligned_path <- file.path(phylo_dir, "Endo_mafft_aligned.fasta")
if (!file.exists(endo_aligned_path)) {
  stop("Endo_mafft_aligned.fasta not found in phylo_dir. This script expects the ",
       "same aligned FASTA used to build Figure 5A's tree (Endo_iqtree.treefile). ",
       "Re-run the alignment step from Endozoicomonas_Phylogenetic_Analysis_FINAL.R ",
       "if this file is missing.")
}
endo_aligned <- read.dna(endo_aligned_path, format = "fasta")
cat("Endozoicomonas alignment loaded:", nrow(endo_aligned), "sequences,",
    ncol(endo_aligned), "aligned positions\n")

# --- Mycoplasma: same alignment used for tree_bs_myco (IQ-TREE input) ---
myco_aligned_path <- file.path(phylo_dir, "Myco_mafft_aligned.fasta")
if (!file.exists(myco_aligned_path)) {
  # Fall back to the R-side msa() alignment saved by the Mycoplasma script,
  # if the IQ-TREE input FASTA isn't found under this name.
  myco_aligned_path <- file.path(myco_output, "Mycoplasma_alignment_CLEAN.fasta")
  cat("Note: Myco_mafft_aligned.fasta not found; falling back to",
      "Mycoplasma_alignment_CLEAN.fasta\n")
}
if (!file.exists(myco_aligned_path)) {
  stop("No Mycoplasma alignment found in phylo_dir or myco_output. Re-run the ",
       "alignment step from Mycoplasma_Phylogenetic_Analysis.R first.")
}
myco_aligned <- read.dna(myco_aligned_path, format = "fasta")
cat("Mycoplasma alignment loaded:", nrow(myco_aligned), "sequences,",
    ncol(myco_aligned), "aligned positions\n")

# =============================================================================
# STEP 2: Compute pairwise base pair (nucleotide) differences
# Raw count of differing sites between each pair of aligned sequences
# (gaps excluded via pairwise.deletion = TRUE, so only positions present
# in both sequences are compared - standard practice for raw bp distance)
# =============================================================================

cat("\n========== STEP 2: Computing pairwise bp differences ==========\n\n")

compute_bp_diff_matrix <- function(aligned_dna) {
  # dist.dna with model="N" returns the raw number of differing sites
  # (not a substitution-corrected distance), i.e. true base pair differences
  d <- dist.dna(aligned_dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
  d
}

endo_bp_matrix <- compute_bp_diff_matrix(endo_aligned)
myco_bp_matrix <- compute_bp_diff_matrix(myco_aligned)

cat("Endozoicomonas bp-difference matrix:", nrow(endo_bp_matrix), "x",
    ncol(endo_bp_matrix), "\n")
cat("Mycoplasma bp-difference matrix:", nrow(myco_bp_matrix), "x",
    ncol(myco_bp_matrix), "\n")

# =============================================================================
# STEP 3: Identify core ASV labels in each alignment
# (reuses the SAME naming convention as tip_data / tip_data_myco in the
# Endo/Myco scripts: "Endo_1|...", "Myco_1|...")
# =============================================================================

cat("\n========== STEP 3: Identifying core ASVs ==========\n\n")

get_core_labels <- function(seq_names, core_prefix, core_ids) {
  core_labels <- seq_names[sapply(core_ids, function(id) {
    grep(paste0("^", core_prefix, "_", id, "\\|"), seq_names)[1]
  })]
  core_labels <- core_labels[!is.na(core_labels)]
  core_labels
}

endo_names <- rownames(endo_aligned)
myco_names <- rownames(myco_aligned)

endo_core <- get_core_labels(endo_names, "Endo", 1:6)
myco_core <- get_core_labels(myco_names, "Myco", 1:10)

cat("Endozoicomonas core ASVs found:", length(endo_core), "/ 6\n")
cat("Mycoplasma core ASVs found:", length(myco_core), "/ 10\n")

# Clean display names: "Endo_1|B. asbestinum+E. caribaeorum" -> "Endo_1"
clean_core_name <- function(x) sub("\\|.*", "", x)

# =============================================================================
# STEP 4: TABLE 1 - Pairwise bp differences among core ASVs themselves
# (the actual request: e.g. Endo_1 vs. Endo_2, Endo_1 vs. Endo_3, etc.,
# and the same for Myco_1 through Myco_10)
# =============================================================================

cat("\n========== STEP 4: Building core-ASV pairwise difference table ==========\n\n")

build_core_pairwise_table <- function(bp_matrix, core_labels, genus_label) {
  n <- length(core_labels)
  rows <- list()
  if (n < 2) return(data.frame())
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      a <- core_labels[i]
      b <- core_labels[j]
      if (a %in% rownames(bp_matrix) && b %in% colnames(bp_matrix)) {
        rows[[length(rows) + 1]] <- data.frame(
          Genus       = genus_label,
          ASV_1       = clean_core_name(a),
          ASV_2       = clean_core_name(b),
          BP_Difference = bp_matrix[a, b],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  bind_rows(rows) %>% arrange(ASV_1, ASV_2)
}

endo_core_pairwise <- build_core_pairwise_table(endo_bp_matrix, endo_core, "Endozoicomonas")
myco_core_pairwise <- build_core_pairwise_table(myco_bp_matrix, myco_core, "Mycoplasma")

table1_full <- bind_rows(endo_core_pairwise, myco_core_pairwise)

cat("Table 1 rows (core ASV x core ASV pairs):", nrow(table1_full), "\n")
print(table1_full)

# Also build a square matrix version (ASV x ASV grid) for each genus - often
# easier to scan visually than a long-format list, especially with only
# 6 (Endo) or 10 (Myco) taxa.
build_square_matrix_df <- function(bp_matrix, core_labels) {
  clean_labels <- clean_core_name(core_labels)
  sub_matrix <- bp_matrix[core_labels, core_labels]
  rownames(sub_matrix) <- clean_labels
  colnames(sub_matrix) <- clean_labels
  as.data.frame(sub_matrix) %>%
    rownames_to_column("ASV") %>%
    mutate(across(-ASV, as.integer))
}

endo_square <- build_square_matrix_df(endo_bp_matrix, endo_core)
myco_square <- build_square_matrix_df(myco_bp_matrix, myco_core)

cat("\nEndozoicomonas core ASV bp-difference matrix:\n")
print(endo_square)
cat("\nMycoplasma core ASV bp-difference matrix:\n")
print(myco_square)

# =============================================================================
# STEP 5: TABLE 2 - Each core ASV vs. ALL other ASVs of the same genus
# in the tree (not just the other core ASVs), summarized as
# min/max/mean/SD bp difference - broader context beyond the core-only
# comparison in Table 1
# =============================================================================

cat("\n========== STEP 5: Building core-ASV vs. all-ASVs summary table ==========\n\n")

build_core_vs_all_summary <- function(bp_matrix, all_names, core_labels, genus_label) {
  # "All ASVs" = every study sequence in this alignment (excludes REF| entries),
  # so this shows how divergent each core ASV is from the full range of other
  # ASVs of the same genus included in the tree - not just the other 5-9 core
  # ASVs, but the complete set.
  study_seqs <- all_names[!grepl("^REF\\|", all_names)]
  
  rows <- list()
  for (core in core_labels) {
    if (!(core %in% rownames(bp_matrix))) next
    other_seqs <- setdiff(study_seqs, core)
    diffs <- bp_matrix[core, other_seqs]
    diffs <- diffs[!is.na(diffs)]
    if (length(diffs) == 0) next
    rows[[length(rows) + 1]] <- data.frame(
      Genus           = genus_label,
      Core_ASV        = clean_core_name(core),
      N_ASVs_Compared = length(diffs),
      Min_BP_Diff     = min(diffs),
      Max_BP_Diff     = max(diffs),
      Mean_BP_Diff    = round(mean(diffs), 1),
      SD_BP_Diff      = round(sd(diffs), 1),
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows) %>% arrange(Core_ASV)
}

endo_all_summary <- build_core_vs_all_summary(endo_bp_matrix, endo_names, endo_core, "Endozoicomonas")
myco_all_summary <- build_core_vs_all_summary(myco_bp_matrix, myco_names, myco_core, "Mycoplasma")

table2_full <- bind_rows(endo_all_summary, myco_all_summary)

cat("Table 2 rows (one per core ASV):", nrow(table2_full), "\n")
print(table2_full)

# =============================================================================
# STEP 6: Save CSVs (for your own records / re-import) and build the
# Word supplementary table document, matching the existing caption and
# table style used in Supplementary_material_Revised_FINAL.docx
# (bold "Table SX." + bold description + table with bold header row,
#  species/taxa names in italics)
# =============================================================================

cat("\n========== STEP 6: Exporting CSVs and Word document ==========\n\n")

write.csv(table1_full, file.path(supp_output, "TableS8a_core_ASV_pairwise_bp_differences.csv"),
          row.names = FALSE)
write.csv(endo_square, file.path(supp_output, "TableS8a_Endo_core_ASV_matrix.csv"),
          row.names = FALSE)
write.csv(myco_square, file.path(supp_output, "TableS8a_Myco_core_ASV_matrix.csv"),
          row.names = FALSE)
write.csv(table2_full, file.path(supp_output, "TableS8b_core_ASV_vs_all_ASV_summary.csv"),
          row.names = FALSE)

# --- Build flextable objects styled to match existing supplementary tables ---
# Table S8a is presented as two square matrices (one per genus) rather than
# a long pairwise list - much easier to scan for 6x6 / 10x10 taxa.

build_matrix_flextable <- function(square_df) {
  square_df %>%
    flextable() %>%
    bold(part = "header") %>%
    italic(j = "ASV", italic = FALSE, part = "body") %>%
    align(align = "center", part = "all") %>%
    align(j = "ASV", align = "left", part = "all") %>%
    bold(j = "ASV", part = "body") %>%
    autofit() %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    border_outer(border = fp_border(color = "black", width = 1)) %>%
    border_inner(border = fp_border(color = "grey70", width = 0.5))
}

ft_table1_endo <- build_matrix_flextable(endo_square)
ft_table1_myco <- build_matrix_flextable(myco_square)

ft_table2 <- table2_full %>%
  flextable() %>%
  set_header_labels(
    Genus = "Genus", Core_ASV = "Core ASV",
    N_ASVs_Compared = "N ASVs compared",
    Min_BP_Diff = "Min", Max_BP_Diff = "Max",
    Mean_BP_Diff = "Mean", SD_BP_Diff = "SD"
  ) %>%
  add_header_row(values = c("", "", "", "BP difference vs. all other ASVs"),
                 colwidths = c(1, 1, 1, 4)) %>%
  italic(j = "Genus", italic = TRUE, part = "body") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Genus", "Core_ASV"), align = "left", part = "all") %>%
  merge_v(j = c("Genus"), part = "body") %>%
  autofit() %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Calibri", part = "all") %>%
  border_outer(border = fp_border(color = "black", width = 1)) %>%
  border_inner_h(border = fp_border(color = "grey70", width = 0.5))

# --- Assemble Word document with captions matching existing table style ---

doc <- read_docx() %>%
  body_add_par("Supplementary material", style = "heading 1") %>%
  body_add_par("") %>%
  
  # Table S8a caption (bold label + bold description, matching e.g. Table S1's style)
  body_add_fpar(fpar(
    ftext("Table S8a. ", fp_text(bold = TRUE)),
    ftext(paste0(
      "Pairwise base pair differences among core Endozoicomonas (Endo_1\u2013Endo_6) ",
      "and Mycoplasma (Myco_1\u2013Myco_10) ASVs."
    ), fp_text(bold = TRUE)),
    ftext(paste0(
      " Values represent raw nucleotide differences (excluding alignment gaps) ",
      "over the aligned V4 region used to construct the phylogenetic trees in ",
      "Figures 5A and 6A."
    ), fp_text(bold = FALSE))
  )) %>%
  body_add_par("") %>%
  body_add_par("Endozoicomonas", style = "Normal") %>%
  body_add_flextable(ft_table1_endo) %>%
  body_add_par("") %>%
  body_add_par("Mycoplasma", style = "Normal") %>%
  body_add_flextable(ft_table1_myco) %>%
  body_add_par("") %>%
  body_add_par("") %>%
  
  # Table S8b caption
  body_add_fpar(fpar(
    ftext("Table S8b. ", fp_text(bold = TRUE)),
    ftext(paste0(
      "Base pair divergence of each core Endozoicomonas and Mycoplasma ASV ",
      "relative to all other ASVs of the same genus included in the ",
      "phylogenetic reconstructions (Figures 5A and 6A)."
    ), fp_text(bold = TRUE)),
    ftext(paste0(
      " Min, Max, Mean, and SD summarize the distribution of pairwise base pair ",
      "differences between the indicated core ASV and every other ASV of the ",
      "same genus in the alignment."
    ), fp_text(bold = FALSE))
  )) %>%
  body_add_par("") %>%
  body_add_flextable(ft_table2)

print(doc, target = file.path(supp_output, "Supplementary_Table_S8_bp_differences.docx"))

cat("\n=============================================================================\n")
cat("DONE\n")
cat("=============================================================================\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n###############################################################\n")
cat("SUPPLEMENTARY MATERIALS COMPLETE\n")
cat("###############################################################\n")
cat("\nFiles saved to:", supp_dir, "\n\n")

cat("Tables:\n")
cat("  Table_S2_sample_counts.docx/.csv\n")
cat("  Table_S3_alpha_diversity.docx/.csv\n")
cat("  Table_S3b_dunn_pairwise.csv\n")
cat("  Table_S4_PERMANOVA.docx/.csv\n")
cat("  Table_S5_pairwise_PERMANOVA.docx/.csv  (includes seawater comparisons)\n")
cat("  Table_S6_ANCOMBC2.docx/.csv\n")
cat("  Table_S7_core_microbiome_all_thresholds.csv\n")
cat("  Table_S7b_bootstrap_validated_core.csv\n")

cat("\nFigures:\n")
cat("  FigS_Species_PCoA_Aitchison_combined.png/.pdf\n")
cat("  FigS_PCoA_Aitchison_[species].png  (4 individual panels)\n")
cat("\nOutput files (in", supp_output, "):\n")
cat("  - TableS8a_core_ASV_pairwise_bp_differences.csv (long-format list)\n")
cat("  - TableS8a_Endo_core_ASV_matrix.csv (6x6 matrix)\n")
cat("  - TableS8a_Myco_core_ASV_matrix.csv (10x10 matrix)\n")
cat("  - TableS8b_core_ASV_vs_all_ASV_summary.csv\n")
cat("  - Supplementary_Table_S8_bp_differences.docx  <-- insert into main supplementary file\n")
cat("\nNext step: copy the two tables from this .docx into ",
    "Supplementary_material__Revised_FINAL.docx, immediately after Table S7,\n",
    "renumbering as needed if additional supplementary tables are added before this point.\n", sep = "")