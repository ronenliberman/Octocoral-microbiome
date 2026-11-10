# =============================================================================
# 02_stable_ASV_IDs.R — One ASV label per ASV, used by every figure and table
# =============================================================================
# Editor comment addressed: "Endo 4, 5 and 6 … you have 1, 2 and 3 twice".
# In R1, Fig 3 numbered ASVs within each host by hash order, Fig 4 by order of
# appearance, and Fig 5 used Endo_1–6 — so "Endozoicomonas 2" meant different
# ASVs in different figures. This script builds ONE lookup table:
#
#   Endozoicomonadaceae ASVs  Endo_1–Endo_6 kept (core, cited in text and S8);
#                             other VERIFIED ASVs numbered Endo_7, Endo_8, …
#                             by total reads.
#   Mycoplasmataceae ASVs     Myco_1–Myco_10 kept (R1 labels, cited in text and S8);
#                             other verified ASVs Myco_11, … by total reads.
#   Unverified ASVs           Unv_1, Unv_2, … by total reads (never called
#                             Endozoicomonas / Mycoplasma in figures).
#   Every other ASV           lowest resolved SILVA rank + a number that is
#                             global across Fig 3 and Fig 4 (by total reads).
#
# Figure_Label is the text to print on figures, e.g.
#   "Endozoicomonas (Endo_2)", "Endozoicomonadaceae (Endo_14)",
#   "Unverified (Unv_3)", "Filomicrobium 2".
#
# Inputs:  out_dir/ASV_validation_status.csv   (from 01)
#          Output/Mycoplasma_phylogeny/Mycoplasma_master_table.csv (R1 Myco labels)
#          Output/Core_Microbiome_AllThresholds_Summary_TSS.csv   (Fig 3 ASVs)
#          Output/ANCOMBC/Figure4_ANCOMBC2_{Year,Site}_significant.csv (Fig 4 ASVs)
# Output:  out_dir/ASV_stable_IDs.csv
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

cat("\n========== 02: Stable ASV IDs ==========\n\n")

d <- load_octocoral_dataset(include_seawater = FALSE)

validation <- read.csv(file.path(out_dir, "ASV_validation_status.csv"), stringsAsFactors = FALSE) %>%
  filter(Query_set %in% c("Endozoicomonadaceae", "Mycoplasmataceae")) %>%
  distinct(ASV_ID, .keep_all = TRUE)

reads <- tibble(ASV_ID = rownames(d$otu_table), Total_Reads = rowSums(d$otu_table))

# =============================================================================
# STEP 1: Endozoicomonadaceae
# =============================================================================

endo <- validation %>%
  filter(Query_set == "Endozoicomonadaceae") %>%
  left_join(reads, by = "ASV_ID") %>%
  mutate(Stable_ID = unname(core_endo_map[ASV_ID]))

endo_rest <- endo %>%
  filter(is.na(Stable_ID), Verified %in% TRUE) %>%
  arrange(desc(Total_Reads)) %>%
  mutate(Stable_ID = paste0("Endo_", row_number() + length(core_endo_map)))

endo <- endo %>%
  filter(!is.na(Stable_ID)) %>%
  bind_rows(endo_rest) %>%
  bind_rows(endo %>% filter(is.na(Stable_ID), !(Verified %in% TRUE)))

# =============================================================================
# STEP 2: Mycoplasmataceae (keep R1 Myco_1–Myco_10)
# =============================================================================

myco_r1 <- read.csv(file.path(blast_dir, "Mycoplasma_master_table.csv"),
                    stringsAsFactors = FALSE) %>%
  dplyr::select(ASV_ID, R1_Label = Label)
myco_keep <- myco_r1 %>% filter(R1_Label %in% paste0("Myco_", 1:10))

myco <- validation %>%
  filter(Query_set == "Mycoplasmataceae") %>%
  left_join(reads, by = "ASV_ID") %>%
  left_join(myco_keep, by = "ASV_ID") %>%
  mutate(Stable_ID = R1_Label) %>%
  dplyr::select(-R1_Label)

myco_rest <- myco %>%
  filter(is.na(Stable_ID), Verified %in% TRUE) %>%
  arrange(desc(Total_Reads)) %>%
  mutate(Stable_ID = paste0("Myco_", row_number() + 10))

myco <- myco %>%
  filter(!is.na(Stable_ID)) %>%
  bind_rows(myco_rest) %>%
  bind_rows(myco %>% filter(is.na(Stable_ID), !(Verified %in% TRUE)))

# =============================================================================
# STEP 3: Unverified ASVs (either family) get Unv_ numbers
# =============================================================================

genus_ids <- bind_rows(endo, myco)
unv <- genus_ids %>%
  filter(is.na(Stable_ID)) %>%
  arrange(desc(Total_Reads)) %>%
  mutate(Stable_ID = paste0("Unv_", row_number()))
genus_ids <- genus_ids %>% filter(!is.na(Stable_ID)) %>% bind_rows(unv)

# ---- Class-level name for unverified ASVs -------------------------------------
# If the best core_nt hit (01b; any hit, incl. MAG genomes) is >= 78.5% identical
# (Yarza et al. 2014 class threshold) over >= 90% of the ASV and its NCBI lineage
# has a class, the ASV is labelled by that class, e.g. "c_Spirochaetia (Unv_1)".

PIDENT_CLASS <- 78.5
nt_top10 <- file.path(blast_dir, "core_nt_hits_top10.csv")
unv_class <- tibble(ASV_ID = character(), Reassigned_class = character(),
                    Class_hit = character(), Class_hit_pident = numeric())
if (file.exists(nt_top10)) {
  best_nt <- read.csv(nt_top10, stringsAsFactors = FALSE) %>%
    filter(ASV_ID %in% unv$ASV_ID, Hit_rank == 1)
  taxids <- unique(na.omit(as.character(best_nt$NCBI_taxid)))
  x <- xml2::read_xml(httr::content(httr::GET(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    query = list(db = "taxonomy", id = paste(taxids, collapse = ","), retmode = "xml")),
    "text", encoding = "UTF-8"))
  cls <- map_dfr(xml2::xml_find_all(x, "/TaxaSet/Taxon"), function(t) {
    ranks <- xml2::xml_text(xml2::xml_find_all(t, "./LineageEx/Taxon/Rank"))
    nms   <- xml2::xml_text(xml2::xml_find_all(t, "./LineageEx/Taxon/ScientificName"))
    tibble(NCBI_taxid = xml2::xml_text(xml2::xml_find_first(t, "./TaxId")),
           NCBI_class = if (any(ranks == "class")) nms[ranks == "class"][1] else NA_character_)
  })
  unv_class <- best_nt %>%
    mutate(NCBI_taxid = as.character(NCBI_taxid)) %>%
    left_join(cls, by = "NCBI_taxid") %>%
    filter(!is.na(NCBI_class), pident >= PIDENT_CLASS, coverage_pct >= 90) %>%
    transmute(ASV_ID, Reassigned_class = NCBI_class,
              Class_hit = sub(" 16S.*| genome.*", "", NCBI_title), Class_hit_pident = pident)
  cat("Unverified ASVs named at class level:", nrow(unv_class), "of", nrow(unv), "\n")
  print(unv_class %>% left_join(unv %>% dplyr::select(ASV_ID, Stable_ID), by = "ASV_ID") %>%
          dplyr::select(Stable_ID, Reassigned_class, Class_hit_pident, Class_hit))
}

genus_ids <- genus_ids %>%
  left_join(unv_class, by = "ASV_ID") %>%
  mutate(
    Figure_Label = case_when(
      Validation_status == "Verified (genus)"  ~ paste0(SILVA_genus, " (", Stable_ID, ")"),
      Validation_status == "Verified (family)" ~ paste0(SILVA_family, " (", Stable_ID, ")"),
      !is.na(Reassigned_class)                 ~ paste0("c_", Reassigned_class, " (", Stable_ID, ")"),
      # Core/R1-labelled ASVs that fail keep their ID but are flagged
      TRUE ~ paste0("Unverified (", Stable_ID, ")")
    )
  ) %>%
  dplyr::select(ASV_ID, Query_set, Stable_ID, Figure_Label, Validation_status,
                NCBI_organism, NCBI_pident, SILVA_family, SILVA_genus,
                Reassigned_class, Class_hit, Class_hit_pident, Total_Reads)

flagged <- genus_ids %>% filter(grepl("^(Endo|Myco)_", Stable_ID), !grepl("^Verified", Validation_status))
if (nrow(flagged) > 0) {
  cat("WARNING: labelled ASVs that did NOT pass validation (decide how to report):\n")
  print(flagged %>% dplyr::select(Stable_ID, NCBI_organism, NCBI_pident, Validation_status))
}

# =============================================================================
# STEP 4: All other ASVs shown in Fig 3 or Fig 4 — global numbering
# =============================================================================

fig3_asvs <- read.csv(file.path(old_out_dir, "Core_Microbiome_AllThresholds_Summary_TSS.csv"),
                      stringsAsFactors = FALSE)$ASV_ID
# Fig 4 ASVs: the PEE2 re-run (06c) when it exists, otherwise the R1 results
pee2_dir   <- file.path(out_dir, "ANCOMBC2_PEE2")
fig4_names <- c("Figure4_ANCOMBC2_Year_significant.csv", "Figure4_ANCOMBC2_Site_significant.csv")
fig4_files <- if (all(file.exists(file.path(pee2_dir, fig4_names)))) file.path(pee2_dir, fig4_names) else
  file.path(old_out_dir, "ANCOMBC", fig4_names)
cat("Fig 4 ASVs taken from:", dirname(fig4_files[1]), "\n")
fig4_asvs <- unlist(lapply(fig4_files, function(f) str_extract(read.csv(f)$Taxa, "[0-9a-f]{32}")))
fig_asvs <- setdiff(unique(c(fig3_asvs, fig4_asvs)), genus_ids$ASV_ID)

lowest_label <- function(tax_row) {
  ranks <- c("Genus", "Family", "Order", "Class", "Phylum")
  prefix <- c(Genus = "", Family = "f_", Order = "o_", Class = "c_", Phylum = "p_")
  for (r in ranks) {
    v <- strip_prefix(tax_row[[r]])
    if (!is.na(v) && v != "" && !grepl("unclassified|uncultured", v, ignore.case = TRUE)) {
      return(if (r == "Genus") v else paste0("Uncl. ", v))
    }
  }
  "Uncl. Bacteria"
}

# Labels corrected after the independent audit (Output/Audit/Audit_Status_Report.md).
# Applied BEFORE numbering, so labels sharing a Display are renumbered consistently.
label_overrides <- tibble::tribble(
  ~ASV_ID,                             ~Display,                ~Evidence,
  "0b5cf529860df00b9544cd19df6c73b9",  "Uncl. Mollicutes",
  "SILVA: Candidatus Hepatoplasma (Entomoplasmatales Incertae Sedis). No Hepatoplasma among 50 core_nt hits; best hit 92.3% (Mycoplasmatales MAG), best genus-level hit 90.7% (Ca. Mariniplasma, Mycoplasmataceae); type strains <= 80.7%; dada2/RDP genus bootstrap 36%. Methods rule supports class only. Audit A2/A3/A7."
)
# Not overridden (checked): 99ad7828 'Uncl. Rhizobiaceae 1'. V4 is 98% to Rhodobacterales type strains, but
# 100% to cultured strain SCSIO 12814 (MZ262827), whose full-length 16S is 93.6-95.5% to Hyphomicrobiales
# type strains only (Salaquimonas, Pseudahrensia, Mesorhizobium) — the SILVA label stands.

tax_df <- d$tax_table %>% as.data.frame()
other_ids <- tibble(ASV_ID = fig_asvs) %>%
  filter(ASV_ID %in% rownames(tax_df)) %>%
  mutate(Display = sapply(ASV_ID, function(a) lowest_label(tax_df[a, ])),
         Display = ifelse(ASV_ID %in% label_overrides$ASV_ID,
                          label_overrides$Display[match(ASV_ID, label_overrides$ASV_ID)], Display)) %>%
  left_join(reads, by = "ASV_ID") %>%
  group_by(Display) %>%
  arrange(desc(Total_Reads), .by_group = TRUE) %>%
  mutate(Figure_Label = if (n() > 1) paste(Display, row_number()) else Display) %>%
  ungroup() %>%
  transmute(ASV_ID, Query_set = "Other", Stable_ID = Figure_Label, Figure_Label,
            Validation_status = NA_character_, Total_Reads)

stable_ids <- bind_rows(genus_ids, other_ids)
write.csv(stable_ids, file.path(out_dir, "ASV_stable_IDs.csv"), row.names = FALSE)

cat("\nStable IDs written:", nrow(stable_ids), "ASVs ->", file.path(out_dir, "ASV_stable_IDs.csv"), "\n")
stable_ids %>% count(Query_set) %>% print()
stopifnot(!any(duplicated(stable_ids$Figure_Label)))
cat("Check passed: no Figure_Label is used for two different ASVs.\n")
