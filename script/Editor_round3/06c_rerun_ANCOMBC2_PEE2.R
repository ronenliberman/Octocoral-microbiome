# =============================================================================
# 06c_rerun_ANCOMBC2_PEE2.R — Re-run ANCOM-BC2 (Figure 4) on the PEE2 table
# =============================================================================
# Decision (Ronen, 1 Oct 2026): re-run ANCOM-BC2 on the same ASV table as all
# other analyses. The R1 Figure 4 was fitted on
#   16 rRNA analysis/raw data/table-NoEMC.qza (+ taxonomy.qza)
# and removed only 2 of the 3 excluded samples (EFLE_DC3_3_2024 was kept).
#
# This runner executes script/Figure4_species_ANCOMBC_REVISION2.R VERBATIM
# (same model, settings and comparisons) with only these substitutions:
#   - input: raw files/table-NoEMC-PEE2.qza, taxonomy-PEE2.qza, metadata_all_edited.csv
#   - samples_to_remove: the standard three samples
#   - output_dir: output/Editor_round3/ANCOMBC2_PEE2/  (R1 files untouched)
#   - the plotting part is skipped here (06b draws Figure 4 from these results
#     after 02 has assigned stable labels to any newly significant ASVs)
#
# Run order afterwards: 02_stable_ASV_IDs.R -> 05 -> 06b -> 08
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

ancom_out <- file.path(out_dir, "ANCOMBC2_PEE2")
dir.create(ancom_out, showWarnings = FALSE, recursive = TRUE)

src <- readLines(file.path(rev_dir, "script", "Figure4_species_ANCOMBC_REVISION2.R"), warn = FALSE)

# keep everything up to (not including) STEP 7 — model fitting + saving results
stop_at <- grep("STEP 7 \\(REVISED\\)", src)[1] - 2
code <- src[1:stop_at]

sub_once <- function(code, pattern, replacement) {
  hit <- grep(pattern, code, fixed = TRUE)
  if (length(hit) != 1) stop("Expected exactly one match for: ", pattern, " (found ", length(hit), ")")
  code[hit] <- sub(pattern, replacement, code[hit], fixed = TRUE)
  code
}

code <- code[!grepl("^setwd\\(", code)]
code <- sub_once(code, '"raw data/table-NoEMC.qza"',
                 sprintf('"%s"', file.path(raw_dir, "table-NoEMC-PEE2.qza")))
code <- sub_once(code, '"raw data/taxonomy.qza"',
                 sprintf('"%s"', file.path(raw_dir, "taxonomy-PEE2.qza")))
code <- sub_once(code, '"raw data/metadata_all_edited.csv"',
                 sprintf('"%s"', file.path(raw_dir, "metadata_all_edited.csv")))
code <- sub_once(code, 'samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024")',
                 'samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024", "EFLE_DC3_3_2024")')
# output_dir is defined once near the top; point it (and the dir.create) at ancom_out
code <- sub_once(code, 'output_dir <- file.path(here::here(), "output")',
                 'output_dir <- ancom_out')

cat("Running ANCOM-BC2 on PEE2 (", length(code), "lines of the R1 script) ...\n")
t0 <- Sys.time()
eval(parse(text = code), envir = globalenv())
cat("\nANCOM-BC2 finished in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# ---- compare with R1 -----------------------------------------------------------
# The R1 result files are not in the repository, so the comparison runs only where
# they are present (it was made on 1 Oct 2026: the same 24 year and 9 site results).
r1_year_csv <- file.path(old_out_dir, "ANCOMBC", "Figure4_ANCOMBC2_Year_significant.csv")
r1_site_csv <- file.path(old_out_dir, "Figure4_ANCOMBC2_Site_significant.csv")
if (file.exists(r1_year_csv) && file.exists(r1_site_csv)) {
old_year <- read.csv(r1_year_csv)
old_site <- read.csv(r1_site_csv)
key <- function(df) paste(df$Species, df$Comparison, stringr::str_extract(df$Taxa, "[0-9a-f]{32}"))
cmp <- function(new, old, what) {
  kn <- key(new); ko <- key(old)
  cat(sprintf("%-5s R1 %2d | PEE2 %2d | shared %2d | new only %2d | lost %2d\n",
              what, length(ko), length(kn), length(intersect(kn, ko)),
              length(setdiff(kn, ko)), length(setdiff(ko, kn))))
}
cat("\nSignificant ASV x comparison (R1 vs PEE2):\n")
cmp(sig_results, old_year, "Year")
cmp(sig_results_site, old_site, "Site")
}
cat("\nResults in:", ancom_out, "\nNext: 02_stable_ASV_IDs.R, then 05, 06b, 08.\n")
