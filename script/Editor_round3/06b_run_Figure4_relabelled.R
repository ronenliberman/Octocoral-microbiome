# =============================================================================
# 06b_run_Figure4_relabelled.R — Redraw Figure 4 from the SAVED ANCOM-BC2 results
# =============================================================================
# ANCOM-BC2 is not re-run. This runner:
#   1. reads the saved significant-ASV tables (same numbers as the R1 figure),
#   2. takes the plotting code (STEP 7 to the end) VERBATIM from
#      script/Figure4_species_ANCOMBC_REVISION2.R,
#   3. inserts apply_stable_labels() (06_Figure4_relabel.R) right after STEP 7C,
#   4. writes every output into Output/Editor_round3/ (R1 files untouched).
# The manuscript figure is the bubble plot:
#   Output/Editor_round3/Figure4_ANCOMBC2_bubble_ALTERNATIVE_revision2.png
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))
source(file.path(script_dir, "06_Figure4_relabel.R"))

suppressPackageStartupMessages({
  library(cowplot)
  library(ggh4x)
  library(stringr)
})

cat("\n========== 06b: Figure 4 with stable ASV labels ==========\n\n")

# --- saved ANCOM-BC2 results -------------------------------------------------
# Site: the 17 Aug 2026 file (site names BrS/BrD, as in the submitted figure);
# identical statistics to Output/ANCOMBC/ (19 Mar), which still says DfS/DfD.
# Preferred: the PEE2 re-run from 06c (decision 1 Oct 2026); fallback: R1 results.
pee2_dir <- file.path(out_dir, "ANCOMBC2_PEE2")
if (all(file.exists(file.path(pee2_dir, c("Figure4_ANCOMBC2_Year_significant.csv",
                                          "Figure4_ANCOMBC2_Site_significant.csv"))))) {
  cat("Using ANCOM-BC2 results re-run on PEE2 (06c)\n")
  sig_results      <- read.csv(file.path(pee2_dir, "Figure4_ANCOMBC2_Year_significant.csv"), stringsAsFactors = FALSE)
  sig_results_site <- read.csv(file.path(pee2_dir, "Figure4_ANCOMBC2_Site_significant.csv"), stringsAsFactors = FALSE)
} else {
  cat("Using R1 ANCOM-BC2 results\n")
  sig_results      <- read.csv(file.path(old_out_dir, "ANCOMBC", "Figure4_ANCOMBC2_Year_significant.csv"),
                               stringsAsFactors = FALSE)
  sig_results_site <- read.csv(file.path(old_out_dir, "Figure4_ANCOMBC2_Site_significant.csv"),
                               stringsAsFactors = FALSE)
}
cat("Significant rows: year", nrow(sig_results), "| site", nrow(sig_results_site), "\n")

# --- plotting code from the R1 script, verbatim --------------------------------
src <- readLines(file.path(rev_dir, "script", "Figure4_species_ANCOMBC_REVISION2.R"), warn = FALSE)
start <- grep("STEP 7 \\(REVISED\\)", src)[1] - 1
stopifnot(!is.na(start))
code <- src[start:length(src)]

# insert the relabelling right before the within-facet de-duplication (end of 7C)
ins <- grep("^# Handle any remaining duplicates within a single facet cell", code)[1]
stopifnot(!is.na(ins))
code <- append(code, c(
  "# ---- round 3: stable ASV labels (06_Figure4_relabel.R) ----",
  "plot_data_year <- apply_stable_labels(plot_data_year)",
  "plot_data_site <- apply_stable_labels(plot_data_site)"
), after = ins - 1)

# panel letters: the R1 code drew "(A" and "(B" (no closing parenthesis)
code <- gsub('draw_label("(A"', 'draw_label("A"', code, fixed = TRUE)
code <- gsub('draw_label("(B"', 'draw_label("B"', code, fixed = TRUE)

# drop the trailing interactive line (table(cut(...))) — not needed for the figure
code <- code[!grepl("^table\\(cut\\(abs\\(bubble_data_year", code)]

output_dir <- out_dir   # every ggsave/write.csv in the R1 code now targets Editor_round3
eval(parse(text = code), envir = globalenv())

cat("\nLabels now used in Figure 4 (Endozoicomonas / Mycoplasma / Unverified rows):\n")
print(bind_rows(plot_data_year %>% mutate(Panel = "A (year)"),
                plot_data_site %>% mutate(Panel = "B (site)")) %>%
        filter(grepl("Endo_|Myco_|Unv_", ASV_label)) %>%
        dplyr::select(Panel, Species, Comparison, ASV_label, lfc) %>%
        arrange(Panel, Species), n = Inf)
