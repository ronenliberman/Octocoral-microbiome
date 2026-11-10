# =============================================================================
# 04_Myco_tree_revised.R — Figure 6 (tree only) for Mycoplasmataceae ASVs
# =============================================================================
# Same treatment as 03 (every tip labelled, larger symbols, no network,
# alignment statistics). Host colours use the 1% rule from R1: a sample
# counts toward an ASV's host only if Mycoplasma is > 1% of its reads.
#
# Also fixes the R1 legend error: the submitted Fig 6 legend listed
# "E. caribaeorum" for the purple symbols, which are E. flexuosa. Here the
# legend is generated from the data with the shared palette (00_setup.R).
#
# IMPORTANT (see Editor_Revision_Fix_Plan.md, section 4.3): with the 16S
# type-strain database, Myco_1 (85.2%), Myco_2 (84.4%) and Myco_7 (84.0%)
# fall below 86.5%. Their best hits are still Mycoplasma, i.e. they look like
# novel lineages. Until the core_nt check (01, Tier B) is done, Myco_1–Myco_10
# are kept in the tree and flagged; set KEEP_LABELLED_MYCO <- FALSE to apply
# the rule strictly.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

KEEP_LABELLED_MYCO <- FALSE   # all of Myco_1–10 pass validation (genus level via core_nt)

cat("\n========== 04: Mycoplasmataceae tree ==========\n\n")

prefix    <- "Myco_r3"
aln_file  <- file.path(mac_dir, paste0(prefix, "_mafft_aligned.fasta"))
tree_file <- file.path(mac_dir, paste0(prefix, "_iqtree.treefile"))
iq_file   <- file.path(mac_dir, paste0(prefix, "_iqtree.iqtree"))

# =============================================================================
# PHASE 1: select tips and write FASTA
# =============================================================================

d        <- load_octocoral_dataset(include_seawater = FALSE)
rep_seqs <- load_rep_seqs()
ids      <- read.csv(file.path(out_dir, "ASV_stable_IDs.csv"), stringsAsFactors = FALSE) %>%
  filter(Query_set == "Mycoplasmataceae")

master <- build_asv_master(d, ids$ASV_ID, group_asvs = ids$ASV_ID,
                           min_group_fraction = MYCO_MIN_FRACTION) %>%
  left_join(ids %>% dplyr::select(-Total_Reads), by = "ASV_ID") %>%
  mutate(Is_Core  = Stable_ID %in% paste0("Myco_", 1:10),
         Verified = grepl("^Verified", Validation_status))

tips <- master %>%
  filter((Verified | (KEEP_LABELLED_MYCO & Is_Core)),
         N_Samples >= MIN_SAMPLES_TREE | Is_Core)

cat("Mycoplasmataceae ASVs in dataset:", nrow(master), "\n")
cat("  verified:", sum(master$Verified), "| unverified:", sum(!master$Verified), "\n")
cat("  labelled Myco_1–10 that fail the rule:",
    paste(master$Stable_ID[master$Is_Core & !master$Verified], collapse = ", "), "\n")
cat("  tree tips:", nrow(tips), "\n")
write.csv(master, file.path(out_dir, "Myco_all_ASVs_with_status.csv"), row.names = FALSE)

# RefSeq NR_ type strains (00_setup.R: myco_references). All four R1 references
# were other bacteria (e.g. the "Acholeplasma" outgroup was Bordetella).
ref <- get_reference_v4(myco_references)
refs <- ref$seqs
ref_map <- ref$table %>% transmute(Display = Label, Safe, Accession, Record_title, V4_length)
cat("Reference sequences:", length(refs), "(", paste(ref_map$Display, collapse = "; "), ")\n")

asv_seqs <- rep_seqs[tips$ASV_ID]
names(asv_seqs) <- tips$Stable_ID
fasta_path <- file.path(mac_dir, paste0(prefix, "_for_tree.fasta"))
writeXStringSet(c(asv_seqs, refs), fasta_path)
write.csv(ref_map, file.path(mac_dir, paste0(prefix, "_reference_names.csv")), row.names = FALSE)

if (file.exists(aln_file) &&
    !setequal(names(readDNAStringSet(aln_file)), names(readDNAStringSet(fasta_path)))) {
  old <- c(aln_file, tree_file, iq_file); file.remove(old[file.exists(old)])
  cat("Old alignment/tree removed (built from a different sequence set).\n")
}

if (!file.exists(tree_file)) {
  cat("\n=== MAC STEP ===\n")
  cat(sprintf('  mafft --localpair --maxiterate 1000 "%s_for_tree.fasta" > "%s_mafft_aligned.fasta"\n', prefix, prefix))
  cat(sprintf('  iqtree3 -s "%s_mafft_aligned.fasta" -m MFP -B 1000 -alrt 1000 -T AUTO --prefix "%s_iqtree"\n', prefix, prefix))
  cat("Copy the alignment, .treefile and .iqtree back into:\n  ", mac_dir, "\nthen run again.\n")
  message("
Phase 1 finished (this is a planned pause, not an error). Run the Mac step, then run this script again.")
  if (interactive()) invokeRestart("abort") else quit(save = "no", status = 0)
}

# =============================================================================
# PHASE 2: alignment statistics + figure
# =============================================================================

aln_stats <- alignment_stats(aln_file, "Mycoplasmataceae") %>%
  bind_cols(iqtree_summary(iq_file))
print(t(aln_stats))
write.csv(aln_stats, file.path(out_dir, "Alignment_stats_Myco.csv"), row.names = FALSE)

tree <- read.tree(tree_file)

tip_data <- tibble(label = tree$tip.label) %>%
  mutate(Is_Ref = grepl("^REF_", label)) %>%
  left_join(tips %>% dplyr::select(Stable_ID, Dominant_Host, Is_Core, Verified),
            by = c("label" = "Stable_ID")) %>%
  left_join(ref_map, by = c("label" = "Safe")) %>%
  mutate(
    # Flag labelled ASVs that are below the family identity threshold
    Tip_Label  = case_when(Is_Ref ~ Display,
                           !Verified ~ paste0(label, "†"),
                           TRUE ~ label),
    Host_Group = ifelse(Is_Ref, "Reference", Dominant_Host),
    Is_Core    = coalesce(Is_Core, FALSE)
  )
stopifnot(!any(is.na(tip_data$Host_Group)))

# Final display (4 Oct 2026): SH-aLRT/UFBoot values at the nodes (see 03 / 00_setup.R)
fig_h <- max(8, 0.30 * nrow(tip_data) + 2)
p <- plot_ml_tree_support_values(tree, tip_data, outgroup_pattern = "Acholeplasma",
                                 fig_w = 11, fig_h = fig_h,
                                 legend_y = 0.95, tag = "Figure 6")   # legend top right
ggsave(file.path(out_dir, "Figure6_Mycoplasmataceae_tree.png"), p,
       width = 11, height = fig_h, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure6_Mycoplasmataceae_tree.pdf"), p,
       width = 11, height = fig_h, bg = "white")

cat("\nFigure 6 saved (", nrow(tip_data), "tips).\n")
cat("† = best type-strain hit below", PIDENT_FAMILY, "% (explain in legend if kept).\n")
