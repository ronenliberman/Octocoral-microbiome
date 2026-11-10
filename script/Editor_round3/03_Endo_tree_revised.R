# =============================================================================
# 03_Endo_tree_revised.R — Figure 5 (tree only) for Endozoicomonadaceae ASVs
# =============================================================================
# Editor comments addressed:
#   - clades missing taxa names  -> every tip is labelled (stable IDs from 02)
#   - host/source symbols too small -> larger points and legend keys
#   - network redundant           -> no panel B; tree only
#   - tree squeezed to the left   -> unverified long-branch ASVs removed,
#                                    tree drawn across the page width
#   - unambiguously aligned bp?   -> alignment statistics written to a table
#
# Two phases (run the script twice):
#   PHASE 1 (no alignment yet): writes the FASTA and prints the Mac commands.
#   PHASE 2 (alignment + IQ-TREE files copied back into mac_dir): alignment
#            statistics and the figure.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

cat("\n========== 03: Endozoicomonadaceae tree ==========\n\n")

prefix   <- "Endo_r3"
aln_file <- file.path(mac_dir, paste0(prefix, "_mafft_aligned.fasta"))
tree_file <- file.path(mac_dir, paste0(prefix, "_iqtree.treefile"))
iq_file  <- file.path(mac_dir, paste0(prefix, "_iqtree.iqtree"))

# =============================================================================
# PHASE 1: select tips and write FASTA
# =============================================================================

d        <- load_octocoral_dataset(include_seawater = FALSE)
rep_seqs <- load_rep_seqs()
ids      <- read.csv(file.path(out_dir, "ASV_stable_IDs.csv"), stringsAsFactors = FALSE) %>%
  filter(Query_set == "Endozoicomonadaceae")

master <- build_asv_master(d, ids$ASV_ID) %>%
  left_join(ids %>% dplyr::select(-Total_Reads), by = "ASV_ID") %>%
  mutate(Is_Core = Stable_ID %in% core_endo_map,
         Verified = grepl("^Verified", Validation_status))

tips <- master %>% filter(Verified, N_Samples >= MIN_SAMPLES_TREE | Is_Core)

cat("Endozoicomonadaceae ASVs in dataset:", nrow(master), "\n")
cat("  verified:", sum(master$Verified), "| unverified:", sum(!master$Verified), "\n")
cat("  tree tips (verified, >=", MIN_SAMPLES_TREE, "samples, or core):", nrow(tips), "\n")
write.csv(master, file.path(out_dir, "Endo_all_ASVs_with_status.csv"), row.names = FALSE)

# Type strains + outgroup: RefSeq NR_ records, downloaded and trimmed to the
# V4 amplicon (00_setup.R: endo_references, get_reference_v4). The R1 GenBank
# references included three sequences that were not Endozoicomonas.
ref <- get_reference_v4(endo_references)
refs <- ref$seqs
ref_map <- ref$table %>% transmute(Display = Label, Safe, Accession, Record_title, V4_length)
cat("Reference sequences:", length(refs), "| V4 lengths:", paste(range(width(refs)), collapse = "-"), "bp\n")

asv_seqs <- rep_seqs[tips$ASV_ID]
names(asv_seqs) <- tips$Stable_ID
fasta_path <- file.path(mac_dir, paste0(prefix, "_for_tree.fasta"))
writeXStringSet(c(asv_seqs, refs), fasta_path)
write.csv(ref_map, file.path(mac_dir, paste0(prefix, "_reference_names.csv")), row.names = FALSE)

# Stale-result guard: remove an old alignment/tree built from a different sequence set
if (file.exists(aln_file) &&
    !setequal(names(readDNAStringSet(aln_file)), names(readDNAStringSet(fasta_path)))) {
  old <- c(aln_file, tree_file, iq_file); file.remove(old[file.exists(old)])
  cat("Old alignment/tree removed (built from a different sequence set).\n")
}

if (!file.exists(tree_file)) {
  cat("\n=== MAC STEP ===\n")
  cat("Copy", paste0(prefix, "_for_tree.fasta"), "to the Mac and run:\n\n")
  cat(sprintf('  mafft --localpair --maxiterate 1000 "%s_for_tree.fasta" > "%s_mafft_aligned.fasta"\n', prefix, prefix))
  cat(sprintf('  iqtree3 -s "%s_mafft_aligned.fasta" -m MFP -B 1000 -alrt 1000 -T AUTO --prefix "%s_iqtree"\n', prefix, prefix))
  cat("\n(mafft --localpair = L-INS-i, as stated in the Methods. Use 'iqtree2' if that is the installed binary.)\n")
  cat("Copy", basename(aln_file), ",", basename(tree_file), "and", basename(iq_file),
      "back into:\n  ", mac_dir, "\nthen run this script again.\n")
  message("
Phase 1 finished (this is a planned pause, not an error). Run the Mac step, then run this script again.")
  if (interactive()) invokeRestart("abort") else quit(save = "no", status = 0)
}

# =============================================================================
# PHASE 2: alignment statistics
# =============================================================================

aln_stats <- alignment_stats(aln_file, "Endozoicomonadaceae") %>%
  bind_cols(iqtree_summary(iq_file))
print(t(aln_stats))
write.csv(aln_stats, file.path(out_dir, "Alignment_stats_Endo.csv"), row.names = FALSE)

# =============================================================================
# PHASE 2: figure
# =============================================================================

tree <- read.tree(tree_file)

tip_data <- tibble(label = tree$tip.label) %>%
  mutate(Is_Ref = grepl("^REF_", label)) %>%
  left_join(tips %>% dplyr::select(Stable_ID, Dominant_Host, Is_Core),
            by = c("label" = "Stable_ID")) %>%
  left_join(ref_map, by = c("label" = "Safe")) %>%
  mutate(
    Tip_Label  = ifelse(Is_Ref, Display, label),
    Host_Group = ifelse(Is_Ref, "Reference", Dominant_Host),
    Is_Core    = coalesce(Is_Core, FALSE)
  )
stopifnot(!any(is.na(tip_data$Host_Group)))

# Final display (4 Oct 2026): SH-aLRT/UFBoot values at the nodes, drawn 16 in tall
# so the dense clades have room (trial: 13_TRIAL_tree_UFBoot_labels.R). The square
# symbols version is plot_ml_tree() alone.
fig_h <- 16
p <- plot_ml_tree_support_values(tree, tip_data, outgroup_pattern = "Kistimonas",
                                 fig_w = 11, fig_h = fig_h, legend_y = 0.5, tag = "Figure 5")
ggsave(file.path(out_dir, "Figure5_Endozoicomonadaceae_tree.png"), p,
       width = 11, height = fig_h, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure5_Endozoicomonadaceae_tree.pdf"), p,
       width = 11, height = fig_h, bg = "white")

cat("\nFigure 5 saved (", nrow(tip_data), "tips, ", fig_h, "in tall).\n")
cat("For the legend/Methods: alignment length", aln_stats$Alignment_length,
    "positions,", aln_stats$Unambiguously_aligned_positions, "unambiguously aligned (gap-free),",
    aln_stats$Parsimony_informative_positions, "parsimony-informative; model", aln_stats$Best_model, "\n")
