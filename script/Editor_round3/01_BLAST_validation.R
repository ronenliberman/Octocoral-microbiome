# =============================================================================
# 01_BLAST_validation.R — Rebuild Appendix 1 and validate Endozoicomonas /
#                          Mycoplasma ASVs by SILVA–NCBI agreement
# =============================================================================
# Editor comment addressed: BLAST hits of ~75% identity make the
# Endozoicomonas classification suspect.
#
# What this script does
#   1. Reads the three raw web-BLAST hit tables already on disk (16S rRNA
#      type-strain database, NR_ accessions):
#        - Output/Endozoicomonas_phylogeny/Endo_blast_results.csv  (all Endozoicomonadaceae ASVs)
#        - Output/Mycoplasma_phylogeny/Myco_blast_results.csv      (all Mycoplasmataceae ASVs)
#        - Output/Editor_round3/ANCOMBC2_r3_16S_hits.csv           (the 26 ASVs with significant
#          ANCOM-BC2 results in the current analysis; written by 01d. Replaces
#          Output/V8941X58016-Alignment-HitTable.csv, the 58 ASVs of an earlier result set)
#   2. Looks up the NCBI lineage of every hit accession BY RANK NAME
#      (fixes the old Appendix 1, where NCBI_Family held the order and
#      NCBI_Genus held the family).
#   3. Applies the validation rule to the best hit of each ASV:
#        Verified (genus)   SILVA family = NCBI family, identity >= 94.5%, SILVA genus = NCBI genus
#        Verified (family)  SILVA family = NCBI family, identity >= 86.5%
#        Unverified         identity < 86.5%, or NCBI family differs from SILVA
#      Thresholds: Yarza et al. 2014 (doi:10.1038/nrmicro3330), set in 00_setup.R.
#   4. Writes a FASTA of Unverified ASVs for a second BLAST against core_nt
#      (to find uncultured relatives missing from the type-strain database),
#      and — if those results exist — applies the same rule to them (Tier B).
#   5. Exports revised Appendix 1 (.csv + .xlsx) with NCBI_ / SILVA_ prefixed
#      column names.
#
# Needs internet (NCBI E-utilities via rentrez). Lineages are cached in
# out_dir/NCBI_lineage_cache.rds so re-runs do not query NCBI again.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

suppressPackageStartupMessages({
  library(rentrez)
  library(XML)
})

# Optional: an NCBI API key raises the rate limit from 3 to 10 requests/s
# set_entrez_key("YOUR_KEY")

cat("\n========== 01: BLAST validation ==========\n\n")

# =============================================================================
# STEP 1: Data needed for SILVA taxonomy and query lengths
# =============================================================================

d        <- load_octocoral_dataset(include_seawater = FALSE)
rep_seqs <- load_rep_seqs()
silva_conf <- load_silva_confidence()

silva_tax <- d$tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  transmute(ASV_ID,
            SILVA_phylum = strip_prefix(Phylum),
            SILVA_class  = strip_prefix(Class),
            SILVA_order  = strip_prefix(Order),
            SILVA_family = strip_prefix(Family),
            SILVA_genus  = strip_prefix(Genus)) %>%
  left_join(silva_conf, by = "ASV_ID")

query_len <- tibble(ASV_ID = names(rep_seqs), Query_length = width(rep_seqs))

# =============================================================================
# STEP 2: Read raw web-BLAST hit tables (no header, 12 standard columns)
# =============================================================================

hit_cols <- c("qseqid", "NCBI_accession", "NCBI_pident", "NCBI_align_length",
              "NCBI_mismatch", "NCBI_gapopen", "NCBI_qstart", "NCBI_qend",
              "NCBI_sstart", "NCBI_send", "NCBI_evalue", "NCBI_bitscore")

read_hit_table <- function(path, set_name) {
  cat("Reading", basename(path), "\n")
  read.csv(path, header = FALSE, col.names = hit_cols, stringsAsFactors = FALSE) %>%
    mutate(
      # qseqid is either the bare hash or "Label|hash|hosts" (truncated at a space)
      ASV_ID = str_extract(qseqid, "[0-9a-f]{32}"),
      Query_set = set_name
    ) %>%
    filter(!is.na(ASV_ID))
}

blast_files <- c(
  Endozoicomonadaceae = file.path(blast_dir, "Endo_blast_results.csv"),
  Mycoplasmataceae    = file.path(blast_dir, "Myco_blast_results.csv"),
  ANCOMBC2            = file.path(blast_dir, "ANCOMBC2_r3_16S_hits.csv")   # from 01d (7 Oct 2026)
)
stopifnot(all(file.exists(blast_files)))

blast_all <- bind_rows(lapply(names(blast_files), function(nm) read_hit_table(blast_files[[nm]], nm)))

# Tier A top-up from 01b: ASVs of the PEE2 table that were not in the R1 BLAST files
for (set in c("Endozoicomonadaceae", "Mycoplasmataceae")) {
  f <- file.path(blast_dir, paste0("TierA_topup_16S_hits_", set, ".csv"))
  if (file.exists(f)) blast_all <- bind_rows(blast_all, read_hit_table(f, set))
}

# Rank hits within each ASV by bitscore, then identity; keep the top 3
blast_top3 <- blast_all %>%
  group_by(Query_set, ASV_ID) %>%
  arrange(desc(NCBI_bitscore), desc(NCBI_pident), .by_group = TRUE) %>%
  distinct(NCBI_accession, .keep_all = TRUE) %>%
  mutate(Hit_rank = row_number()) %>%
  filter(Hit_rank <= 3) %>%
  ungroup() %>%
  left_join(query_len, by = "ASV_ID") %>%
  mutate(NCBI_query_coverage_pct = round(100 * (NCBI_qend - NCBI_qstart + 1) / Query_length, 1))

cat("ASVs with hits:", n_distinct(blast_top3$ASV_ID),
    "| unique accessions:", n_distinct(blast_top3$NCBI_accession), "\n")

# =============================================================================
# STEP 3: NCBI lineage by accession, ranks looked up BY NAME
# =============================================================================
# Adapted from fetch_taxonomy_from_accession() in
# script/NCBI_BLAST_unclassified_ASVs_localBLAST.R (L240–283).

fetch_lineage <- function(accession) {
  empty <- tibble(NCBI_accession = accession, NCBI_taxid = NA_character_,
                  NCBI_organism = NA_character_, NCBI_phylum = NA_character_,
                  NCBI_class = NA_character_, NCBI_order = NA_character_,
                  NCBI_family = NA_character_, NCBI_genus = NA_character_)
  tryCatch({
    taxid <- entrez_summary(db = "nucleotide", id = accession)$taxid
    if (is.null(taxid) || is.na(taxid)) return(empty)
    tax_xml <- xmlParse(entrez_fetch(db = "taxonomy", id = taxid, rettype = "xml"),
                        asText = TRUE)
    lineage <- list()
    for (node in getNodeSet(tax_xml, "//LineageEx/Taxon")) {
      rank <- xmlValue(getNodeSet(node, ".//Rank")[[1]])
      lineage[[rank]] <- xmlValue(getNodeSet(node, ".//ScientificName")[[1]])
    }
    # The organism itself can be the genus (e.g. "Endozoicomonas sp.") — its own rank
    own_rank <- xmlValue(getNodeSet(tax_xml, "/TaxaSet/Taxon/Rank")[[1]])
    own_name <- xmlValue(getNodeSet(tax_xml, "/TaxaSet/Taxon/ScientificName")[[1]])
    lineage[[own_rank]] <- lineage[[own_rank]] %||% own_name
    free(tax_xml)
    tibble(NCBI_accession = accession, NCBI_taxid = as.character(taxid),
           NCBI_organism = own_name,
           NCBI_phylum = lineage[["phylum"]] %||% NA_character_,
           NCBI_class  = lineage[["class"]]  %||% NA_character_,
           NCBI_order  = lineage[["order"]]  %||% NA_character_,
           NCBI_family = lineage[["family"]] %||% NA_character_,
           NCBI_genus  = lineage[["genus"]]  %||% NA_character_)
  }, error = function(e) {
    cat("  FAILED:", accession, conditionMessage(e), "\n")
    empty
  })
}

cache_file <- file.path(blast_dir, "NCBI_lineage_cache.rds")
lineage_cache <- if (file.exists(cache_file)) readRDS(cache_file) else tibble()

get_lineages <- function(accessions) {
  todo <- setdiff(unique(accessions), lineage_cache$NCBI_accession)
  cat("Lineages cached:", nrow(lineage_cache), "| to fetch:", length(todo), "\n")
  for (i in seq_along(todo)) {
    if (i %% 25 == 0) cat("  ", i, "/", length(todo), "\n")
    lineage_cache <<- bind_rows(lineage_cache, fetch_lineage(todo[i]))
    Sys.sleep(0.4)
    if (i %% 50 == 0) saveRDS(lineage_cache, cache_file)
  }
  saveRDS(lineage_cache, cache_file)
  lineage_cache %>% filter(NCBI_accession %in% accessions)
}

lineages <- get_lineages(blast_top3$NCBI_accession)
blast_top3 <- blast_top3 %>% left_join(lineages, by = "NCBI_accession")

# Sanity check the rank fix: Endo_2 (2f40b4ea…) must show Endozoicomonadaceae / Endozoicomonas
chk <- blast_top3 %>% filter(ASV_ID == "2f40b4eaa74d3e96797cbb513a5f43bd", Hit_rank == 1)
cat("Rank check, Endo_2 top hit:", chk$NCBI_family[1], "/", chk$NCBI_genus[1], "\n")

# =============================================================================
# STEP 4: Validation rule (applied to the best hit of each ASV)
# =============================================================================

family_matches <- function(silva_family, ncbi_family) {
  mapply(function(s, n) {
    if (is.na(s) || is.na(n) || s == "") return(FALSE)
    allowed <- family_equivalents[[s]] %||% s
    n %in% allowed
  }, silva_family, ncbi_family, USE.NAMES = FALSE)
}

validate_hits <- function(best_hits) {
  best_hits %>%
    mutate(
      Family_concordant = family_matches(SILVA_family, NCBI_family),
      Genus_concordant  = !is.na(NCBI_genus) & NCBI_genus == SILVA_genus,
      Pident_ge_family  = NCBI_pident >= PIDENT_FAMILY,
      Pident_ge_genus   = NCBI_pident >= PIDENT_GENUS,
      Validation_status = case_when(
        is.na(NCBI_pident)                                          ~ "No BLAST hit",
        Family_concordant & Pident_ge_genus & Genus_concordant      ~ "Verified (genus)",
        Family_concordant & Pident_ge_family                        ~ "Verified (family)",
        !Pident_ge_family                                           ~ "Unverified (identity < 86.5%)",
        TRUE                                                        ~ "Unverified (family discordant)"
      ),
      Verified = grepl("^Verified", Validation_status),
      # Name to use in figures and text
      Reported_taxon = case_when(
        Validation_status == "Verified (genus)"  ~ SILVA_genus,
        Validation_status == "Verified (family)" ~ paste0(SILVA_family, " (unclassified genus)"),
        TRUE ~ paste0("Unverified (SILVA: ", SILVA_genus, ")")
      )
    )
}

best_hits <- blast_top3 %>%
  filter(Hit_rank == 1) %>%
  left_join(silva_tax, by = "ASV_ID")

# ASVs of the two genera that had no BLAST hit at all still need a row
genus_asvs <- silva_tax %>%
  filter(SILVA_family %in% c("Endozoicomonadaceae", "Mycoplasmataceae") |
           grepl("Endozoicomonas|Mycoplasma", SILVA_genus)) %>%
  mutate(Query_set = ifelse(grepl("Mycoplasma", paste(SILVA_family, SILVA_genus)),
                            "Mycoplasmataceae", "Endozoicomonadaceae"))
missing <- genus_asvs %>%
  filter(ASV_ID %in% rownames(d$otu_table),
         !paste(Query_set, ASV_ID) %in% paste(best_hits$Query_set, best_hits$ASV_ID))
if (nrow(missing) > 0) {
  cat("NOTE:", nrow(missing), "Endo/Myco ASVs in the PEE2 table were not in the old BLAST",
      "files (the R1 Endo master used table-NoEMC.qza). They are exported for BLAST below.\n")
}

validation <- validate_hits(best_hits) %>%
  bind_rows(missing %>% mutate(Validation_status = "Not yet BLASTed", Verified = NA)) %>%
  mutate(Validation_tier = "A: 16S type-strain database")

cat("\n=== Validation (Tier A: 16S rRNA type-strain database) ===\n")
validation %>% count(Query_set, Validation_status) %>% print(n = Inf)

# =============================================================================
# STEP 5: Tier B — second BLAST of Unverified ASVs against core_nt
# =============================================================================
# MANUAL STEP (web BLAST, ~5 min):
#   1. Upload out_dir/Unverified_ASVs_for_nt_BLAST.fasta at https://blast.ncbi.nlm.nih.gov/
#      Nucleotide BLAST → Database: "Core nucleotide database (core_nt)",
#      Program: megablast, Max target sequences: 100.
#   2. Download → "Hit Table (CSV)" → save as out_dir/Unverified_ASVs_nt_hits.csv
#   3. Re-run this script. Hits to "uncultured bacterium" (no family in NCBI
#      taxonomy) are skipped; the best hit WITH a family-level lineage is used.

needs_nt <- validation %>%
  filter(Query_set != "ANCOMBC2",
         grepl("^Unverified|Not yet BLASTed", Validation_status))

nt_fasta <- rep_seqs[names(rep_seqs) %in% needs_nt$ASV_ID]
writeXStringSet(nt_fasta, file.path(out_dir, "Unverified_ASVs_for_nt_BLAST.fasta"))
cat("\nFASTA for core_nt BLAST:", length(nt_fasta), "ASVs ->",
    file.path(out_dir, "Unverified_ASVs_for_nt_BLAST.fasta"), "\n")

nt_file <- file.path(blast_dir, "Unverified_ASVs_nt_hits.csv")
if (file.exists(nt_file)) {
  cat("core_nt results found — applying Tier B.\n")
  # Only ASVs that are still unverified after Tier A (the nt file can contain
  # ASVs that the Tier A top-up has since verified)
  nt_hits <- read_hit_table(nt_file, "core_nt") %>% filter(ASV_ID %in% needs_nt$ASV_ID)
  nt_lin  <- get_lineages(nt_hits$NCBI_accession)
  nt_best <- nt_hits %>%
    left_join(nt_lin, by = "NCBI_accession") %>%
    filter(!is.na(NCBI_family)) %>%
    group_by(ASV_ID) %>%
    arrange(desc(NCBI_bitscore), desc(NCBI_pident), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    left_join(silva_tax, by = "ASV_ID") %>%
    left_join(query_len, by = "ASV_ID") %>%
    mutate(NCBI_query_coverage_pct = round(100 * (NCBI_qend - NCBI_qstart + 1) / Query_length, 1),
           Query_set = needs_nt$Query_set[match(ASV_ID, needs_nt$ASV_ID)])
  tier_b <- validate_hits(nt_best) %>%
    mutate(Validation_tier = "B: core_nt (best hit with family-level taxonomy)")
  # Replace Tier A rows only where Tier B verifies the ASV
  rescued <- tier_b %>% filter(Verified)
  cat("Rescued by core_nt:", nrow(rescued), "of", nrow(needs_nt), "\n")
  # Replace Tier A rows by (ASV_ID, Query_set) pair, not by each column separately
  validation <- validation %>%
    anti_join(rescued %>% dplyr::select(ASV_ID, Query_set), by = c("ASV_ID", "Query_set")) %>%
    bind_rows(rescued)
  stopifnot(!any(is.na(validation$Query_set)))
} else {
  cat("No core_nt results yet (", basename(nt_file), ") — Tier A only.\n")
}

# =============================================================================
# STEP 5b: Genus check for family-level ASVs (hits from 01c_core_nt_genus_upgrade.R)
# =============================================================================
# Upgrade "Verified (family)" to "Verified (genus)" when the best core_nt hit
# with genus-level taxonomy has the same family and genus as SILVA at >= 94.5%.
# Never downgrades.

genus_file <- file.path(blast_dir, "Family_level_ASVs_nt_hits.csv")
if (file.exists(genus_file)) {
  fam_level <- validation %>%
    filter(Query_set != "ANCOMBC2", Validation_status == "Verified (family)")
  g_hits <- read_hit_table(genus_file, "core_nt genus check") %>%
    filter(ASV_ID %in% fam_level$ASV_ID)
  g_lin <- get_lineages(g_hits$NCBI_accession)
  g_best <- g_hits %>%
    left_join(g_lin, by = "NCBI_accession") %>%
    filter(!is.na(NCBI_genus)) %>%
    group_by(ASV_ID) %>%
    arrange(desc(NCBI_bitscore), desc(NCBI_pident), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    left_join(silva_tax, by = "ASV_ID") %>%
    left_join(query_len, by = "ASV_ID") %>%
    mutate(NCBI_query_coverage_pct = round(100 * (NCBI_qend - NCBI_qstart + 1) / Query_length, 1),
           Query_set = fam_level$Query_set[match(ASV_ID, fam_level$ASV_ID)])
  upgraded <- validate_hits(g_best) %>%
    filter(Validation_status == "Verified (genus)") %>%
    mutate(Validation_tier = "B: core_nt genus check")
  cat("Upgraded to genus by core_nt:", nrow(upgraded), "of", nrow(fam_level), "family-level ASVs\n")
  validation <- validation %>%
    anti_join(upgraded %>% dplyr::select(ASV_ID, Query_set), by = c("ASV_ID", "Query_set")) %>%
    bind_rows(upgraded)
} else {
  cat("No genus-check results yet (", basename(genus_file), ").\n")
}

# =============================================================================
# STEP 6: Exports
# =============================================================================

appendix_cols <- c(
  "ASV_ID", "Hit_rank",
  "NCBI_accession", "NCBI_organism", "NCBI_pident", "NCBI_align_length",
  "NCBI_query_coverage_pct", "NCBI_mismatch", "NCBI_gapopen", "NCBI_evalue",
  "NCBI_bitscore", "NCBI_phylum", "NCBI_class", "NCBI_order", "NCBI_family",
  "NCBI_genus",
  "SILVA_phylum", "SILVA_class", "SILVA_order", "SILVA_family", "SILVA_genus",
  "SILVA_confidence"
)

make_appendix <- function(set_name) {
  blast_top3 %>%
    filter(Query_set == set_name) %>%
    left_join(silva_tax, by = "ASV_ID") %>%
    left_join(validation %>% filter(Query_set == set_name) %>%
                dplyr::select(ASV_ID, Validation_status, Validation_tier, Reported_taxon),
              by = "ASV_ID") %>%
    dplyr::select(all_of(appendix_cols), Validation_status, Validation_tier, Reported_taxon) %>%
    arrange(ASV_ID, Hit_rank)
}

appendix1 <- list(
  README = tibble(Note = c(
    "Appendix 1. NCBI BLASTn results (top three hits per ASV; web BLAST, 16S rRNA type-strain database).",
    "Columns prefixed NCBI_ come from BLAST and the NCBI Taxonomy lineage of the hit accession; SILVA_ columns are the SILVA 138 classification of the same ASV.",
    "NCBI_pident = % identity over NCBI_align_length; NCBI_query_coverage_pct = aligned query length / ASV length.",
    paste0("Validation (best hit): Verified (genus) = same family and genus, identity >= ", PIDENT_GENUS,
           "%; Verified (family) = same family, identity >= ", PIDENT_FAMILY,
           "%; otherwise Unverified. Thresholds follow Yarza et al. (2014)."),
    "Validation_tier B = verified by a second BLAST against core_nt (uncultured relatives).",
    "Sheet A: all Endozoicomonadaceae and Mycoplasmataceae ASVs used for the phylogenies. Sheet B: the 26 ASVs with significant ANCOM-BC2 results (Figure 4, Table S7); web BLAST of October 2026 against the same database (megablast; blastn for queries without a megablast hit)."
  )),
  A_Endozoicomonadaceae = make_appendix("Endozoicomonadaceae"),
  A_Mycoplasmataceae    = make_appendix("Mycoplasmataceae"),
  B_ANCOMBC2_ASVs       = make_appendix("ANCOMBC2")
)

write_xlsx(appendix1, file.path(out_dir, "Appendix1_BLAST_top3_NCBI_SILVA_revised.xlsx"))
for (nm in names(appendix1)[-1]) {
  write.csv(appendix1[[nm]], file.path(out_dir, paste0("Appendix1_", nm, ".csv")), row.names = FALSE)
}

write.csv(validation, file.path(out_dir, "ASV_validation_status.csv"), row.names = FALSE)

cat("\n=== Summary for the response letter ===\n")
validation %>%
  filter(Query_set != "ANCOMBC2") %>%
  group_by(Query_set) %>%
  summarise(N_ASVs = n(),
            Verified_genus  = sum(Validation_status == "Verified (genus)"),
            Verified_family = sum(Validation_status == "Verified (family)"),
            Unverified      = sum(grepl("^Unverified", Validation_status)),
            No_hit          = sum(Validation_status %in% c("No BLAST hit", "Not yet BLASTed")),
            .groups = "drop") %>%
  print()

cat("\nCore ASVs:\n")
validation %>%
  filter(ASV_ID %in% names(core_endo_map)) %>%
  mutate(Label = core_endo_map[ASV_ID]) %>%
  dplyr::select(Label, NCBI_organism, NCBI_pident, Validation_status) %>%
  arrange(Label) %>% print()

cat("\nFiles written to", out_dir, ":\n",
    " Appendix1_BLAST_top3_NCBI_SILVA_revised.xlsx (+ .csv per sheet)\n",
    " ASV_validation_status.csv\n",
    " Unverified_ASVs_for_nt_BLAST.fasta\n")
