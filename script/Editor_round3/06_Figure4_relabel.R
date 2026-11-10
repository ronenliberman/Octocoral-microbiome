# =============================================================================
# 06_Figure4_relabel.R — Figure 4 (ANCOM-BC2) with the stable ASV labels
# =============================================================================
# ANCOM-BC2 is NOT re-run. Figure 4 in R1 numbered duplicate genera by order
# of appearance within each species ("Endozoicomonas (1)…(5)"), which clashed
# with Endo_1–6 in Fig 5. This file provides apply_stable_labels(), which
# replaces those labels with the Figure_Label from ASV_stable_IDs.csv.
#
# HOW TO USE
#   In script/Figure4_species_ANCOMBC_REVISION2.R, right AFTER STEP 7C
#   (after the two left_join(global_numbering, …) blocks, ~line 488) and
#   BEFORE deduplicate_within_facet(), add:
#
#     source(here::here("script", "Editor_round3", "06_Figure4_relabel.R"))
#     plot_data_year <- apply_stable_labels(plot_data_year)
#     plot_data_site <- apply_stable_labels(plot_data_site)
#
#   then change the ggsave() file names at the end of that script to
#   "Figure4_ANCOMBC2_bubble_round3.png/.pdf" and run it from STEP 7 on
#   (the ANCOM-BC2 results are read back from the saved *_significant.csv
#   files if the model objects are not in memory — see the note below).
#
#   Also change the legend/figure text "Species-specific" -> "Host-specific".
# =============================================================================

.stable_ids_path <- here::here("output", "Editor_round3", "ASV_stable_IDs.csv")

apply_stable_labels <- function(df) {
  ids <- read.csv(.stable_ids_path, stringsAsFactors = FALSE)
  df$ASV_ID_hash <- stringr::str_extract(df$Taxa, "[0-9a-f]{32}")
  lab <- ids$Figure_Label[match(df$ASV_ID_hash, ids$ASV_ID)]
  n_missing <- sum(is.na(lab))
  if (n_missing > 0) {
    message(n_missing, " ASVs not found in ASV_stable_IDs.csv — keeping their R1 label.")
  }
  df$ASV_label <- ifelse(is.na(lab), df$ASV_label, lab)
  df
}

# NOTE — if the ANCOM-BC2 model objects are no longer in memory, rebuild the
# two plotting tables from the saved results instead of re-running the model:
#
#   sig_results      <- read.csv(file.path(output_dir, "ANCOMBC", "Figure4_ANCOMBC2_Year_significant.csv"))
#   sig_results_site <- read.csv(file.path(output_dir, "ANCOMBC", "Figure4_ANCOMBC2_Site_significant.csv"))
#
# and continue from STEP 7A of Figure4_species_ANCOMBC_REVISION2.R.
