# =============================================================================
# 08_supplement_export.R — One supplementary data workbook the editor can open
# =============================================================================
# Editor comment addressed: "the supplemental tables were not included in
# your submission". ScholarOne leaves files designated "Supplementary
# Material" out of the reviewer PDF, and Appendix 1 / Appendix 2 were never
# uploaded at all. This script:
#   1. Computes Table S10: bootstrap frequency of the 90% core ASVs (cited
#      in Methods P36 but never produced as a table).
#   2. Writes ONE workbook, Supplementary_Data_Liberman_et_al.xlsx, with a
#      README sheet, Tables S8–S10, Appendix 1 and Appendix 2.
# Tables S1–S7 stay in the supplementary Word document.
#
# Upload the workbook AND a PDF of the supplementary Word document as
# "Supplementary Material for Review" in ScholarOne, and say so in the
# cover letter.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

cat("\n========== 08: Supplementary workbook ==========\n\n")

# =============================================================================
# STEP 1: Table S10 — bootstrap validation of the 90% core
# =============================================================================
# Decision (Ronen, 1 Oct 2026): sampling WITHOUT replacement (R1 code); the
# Methods are corrected to say "without replacement". Report the result of this
# documented run (E. caribaeorum 22/29) and name the borderline ASVs.
# (A test with replacement also gave 22/29.)

N_BOOT              <- 1000
N_PER_BOOT          <- 20
BOOT_REPLACE        <- FALSE
BOOT_PREV           <- 90
BOOT_ROBUST_PCT     <- 80
BOOT_BORDERLINE_PCT <- 70   # flagged as borderline (not robust) when 70 <= frequency < 80

d <- load_octocoral_dataset(include_seawater = FALSE)
mt_tss <- trans_norm$new(dataset = d)$norm(method = "TSS")
ids <- read.csv(file.path(out_dir, "ASV_stable_IDs.csv"), stringsAsFactors = FALSE)

core_at <- function(otu, prev) {
  p <- rowSums(otu > 0) / ncol(otu) * 100
  names(p)[p >= prev]
}

boot_rows <- list()
for (sp in species_order_levels) {
  otu_sp <- as.matrix(mt_tss$otu_table[, rownames(mt_tss$sample_table)[mt_tss$sample_table$SPP == sp], drop = FALSE])
  otu_sp <- otu_sp[rowSums(otu_sp) > 0, , drop = FALSE]
  full_core <- core_at(otu_sp, BOOT_PREV)
  n_draw <- min(N_PER_BOOT, ncol(otu_sp))
  cat(sp, ":", ncol(otu_sp), "samples,", length(full_core), "ASVs at", BOOT_PREV, "% | drawing", n_draw, "\n")

  set.seed(123)
  hits <- setNames(integer(length(full_core)), full_core)
  for (b in seq_len(N_BOOT)) {
    cols <- sample(colnames(otu_sp), n_draw, replace = BOOT_REPLACE)
    core_b <- core_at(otu_sp[, cols, drop = FALSE], BOOT_PREV)
    hits <- hits + (names(hits) %in% core_b)
  }
  boot_rows[[sp]] <- tibble(
    Species = sp, ASV_ID = full_core,
    Prevalence_pct = round(rowSums(otu_sp[full_core, , drop = FALSE] > 0) / ncol(otu_sp) * 100, 1),
    Bootstrap_frequency_pct = round(100 * hits / N_BOOT, 1)
  ) %>% mutate(Robust_core = Bootstrap_frequency_pct >= BOOT_ROBUST_PCT,
               Borderline  = !Robust_core & Bootstrap_frequency_pct >= BOOT_BORDERLINE_PCT)
}

tableS10 <- bind_rows(boot_rows) %>%
  left_join(ids %>% dplyr::select(ASV_ID, Figure_Label), by = "ASV_ID") %>%
  relocate(Figure_Label, .after = ASV_ID) %>%
  arrange(factor(Species, levels = species_order_levels), desc(Bootstrap_frequency_pct))

cat("\nRobust core (text to report: E. car 22/29, B. asb 6/7, M. mur 2/2, E. flex 1/1):\n")
print(tableS10 %>% group_by(Species) %>%
        summarise(core_90 = n(), robust = sum(Robust_core), borderline = sum(Borderline)))
cat("\nBorderline ASVs (named in Table S10):\n")
print(tableS10 %>% filter(Borderline) %>%
        dplyr::select(Species, Figure_Label, Prevalence_pct, Bootstrap_frequency_pct), n = Inf)

write.csv(tableS10, file.path(out_dir, "TableS10_core90_bootstrap.csv"), row.names = FALSE)

# =============================================================================
# STEP 2: Assemble the workbook
# =============================================================================

read_sheet <- function(f) read.csv(file.path(out_dir, f), stringsAsFactors = FALSE, check.names = FALSE)

sheets_s8 <- lapply(c("S8_pairs", "S8_Endo_matrix", "S8_Myco_matrix", "S9_all_ASVs"),
                    function(s) readxl::read_excel(file.path(out_dir, "TablesS8_S9_revised.xlsx"), sheet = s))
names(sheets_s8) <- c("Table S8 pairs", "Table S8 Endo matrix", "Table S8 Myco matrix", "Table S9")

# Type-strain reference sequences used in the trees (cited in Methods: "accession numbers in Table S9")
ref_tbl <- bind_rows(
  read.csv(file.path(mac_dir, "Endo_r3_reference_names.csv"), stringsAsFactors = FALSE) %>%
    mutate(Tree = "Figure 5 (Endozoicomonadaceae)"),
  read.csv(file.path(mac_dir, "Myco_r3_reference_names.csv"), stringsAsFactors = FALSE) %>%
    mutate(Tree = "Figure 6 (Mycoplasmataceae)")
) %>%
  transmute(Tree, Reference = Display, Accession, NCBI_record_title = Record_title,
            V4_length_bp = V4_length,
            Role = ifelse(grepl("Kistimonas|Acholeplasma", Display), "Outgroup", "Type strain"))
sheets_s8[["Table S9 references"]] <- ref_tbl

app1 <- lapply(c("README", "A_Endozoicomonadaceae", "A_Mycoplasmataceae", "B_ANCOMBC2_ASVs"),
               function(s) readxl::read_excel(file.path(out_dir, "Appendix1_BLAST_top3_NCBI_SILVA_revised.xlsx"), sheet = s))
# Appendix 1 README (written by 01 on 30 Sep) still names SILVA 138.2; the classifier
# was trained on SILVA 138 (qza provenance). 01 is fixed for future runs; patch the copy here.
app1[[1]]$Note <- gsub("SILVA 138.2", "SILVA 138", app1[[1]]$Note, fixed = TRUE)

# 9 Oct 2026: supplementary tables renumbered in order of citation (old S10/S9/S8 -> S3/S4/S5;
# old S3/S5/S4/S6/S7 -> S6/S7/S8/S9/S10). Sheet names follow the new numbers.
names(sheets_s8) <- sub("^Table S9", "Table S4", sub("^Table S8", "Table S5", names(sheets_s8)))
app1[[1]]$Note <- gsub("Table S7", "Table S10", app1[[1]]$Note, fixed = TRUE)

readme <- tibble(Sheet = c(
  "Table S3",
  "Table S4",
  "Table S5 pairs / matrices",
  "Appendix 1 (README, A-Endo, A-Myco, B-ANCOM)",
  "Appendix 2"
), Content = c(
  paste0("Bootstrap validation of the ", BOOT_PREV, "% core (", N_BOOT, " iterations, ", N_PER_BOOT,
         " colonies per iteration, drawn without replacement; seed 123); robust core = present in >= ",
         BOOT_ROBUST_PCT, "% of iterations. Borderline = TRUE marks ASVs present in ", BOOT_BORDERLINE_PCT,
         "-", BOOT_ROBUST_PCT, "% of iterations (not counted as robust; with 1000 iterations the ",
         "Monte Carlo SE is about 1.3 percentage points, so these could fall on either side of the cut-off)."),
  "All Endozoicomonadaceae and Mycoplasmataceae ASVs: occurrence, host, SILVA and NCBI classification, validation status, and whether the ASV is shown in Fig. 5/6. N_Samples = colonies in which the ASV was detected; Mean_RelAbund_Pct = mean relative abundance (%) across those colonies only (not across all colonies, as in Fig. 3). Sheet 'Table S4 references': NCBI RefSeq type-strain sequences (NR_ accessions) used in the phylogenies, trimmed to the 515F-806R region.",
  "Pairwise comparison of core Endozoicomonas (Endo_1–6) and Mycoplasma (Myco_1–10) ASVs over the tree alignments: positions compared, base differences, % identity. Matrices: lower triangle % identity, upper triangle base differences.",
  "NCBI BLASTn top three hits per ASV (16S rRNA type-strain database) with NCBI lineage and SILVA 138 classification; validation status per ASV.",
  "All core ASVs per host at 80, 90 and 100% prevalence; Shown_in_Fig3 marks ASVs with mean relative abundance >= 1%."
))

workbook <- c(
  list(README = readme),
  list("Table S3" = tableS10),
  sheets_s8[c("Table S4", "Table S4 references", "Table S5 pairs", "Table S5 Endo matrix", "Table S5 Myco matrix")],
  setNames(app1, c("Appendix 1 README", "Appendix 1A Endo", "Appendix 1A Myco", "Appendix 1B ANCOM")),
  list("Appendix 2" = read_sheet("Appendix2_core_microbiome_all_thresholds.csv"))
)

out_xlsx <- file.path(out_dir, "Supplementary_Data_Liberman_et_al.xlsx")
write_xlsx(workbook, out_xlsx)
cat("\nWorkbook written:", out_xlsx, "\nSheets:", paste(names(workbook), collapse = " | "), "\n")
