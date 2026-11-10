# =============================================================================
# 01b_core_nt_BLAST_Myco.R — Tier B check for Mycoplasmataceae ASVs
# =============================================================================
# Mycoplasma ASVs that failed the 16S type-strain check (best hit < 86.5%,
# no hit, or best hit outside the family) are BLASTed against NCBI core_nt
# (megablast) through the NCBI BLAST URL API, so uncultured relatives
# (e.g. from corals and other invertebrates) are included.
#
# Rule (same as 01): an ASV is rescued if its best core_nt hit WITH a
# family-level NCBI lineage belongs to Mycoplasmataceae (or a family NCBI
# split from it) at >= 86.5% identity. Hits named "uncultured bacterium"
# carry no family and are skipped for the rule, but still reported.
#
# Needs only the R1 table Mycoplasma_master_BLAST_annotated.csv (sequences
# and type-strain results) and internet access. No QIIME2 files.
# Outputs (Output/Editor_round3/):
#   Myco_core_nt_raw.json                raw BLAST result
#   Myco_core_nt_hits_top10.csv          top 10 hits per ASV with NCBI lineage
#   Myco_core_nt_validation.csv          one row per ASV: Tier A vs Tier B
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(xml2)
})

rev_dir   <- here::here()
out_dir   <- file.path(rev_dir, "output", "Editor_round3")
blast_dir <- file.path(rev_dir, "blast")   # NCBI results are saved here and reused
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PIDENT_FAMILY <- 86.5
PIDENT_GENUS  <- 94.5
myco_families <- c("Mycoplasmataceae", "Metamycoplasmataceae", "Mycoplasmoidaceae")

cat("\n========== 01b: core_nt BLAST of failing Mycoplasma ASVs ==========\n\n")

# =============================================================================
# STEP 1: Queries = ASVs failing Tier A (16S type-strain database)
# =============================================================================

myco <- read.csv(file.path(blast_dir, "Mycoplasma_master_BLAST_annotated.csv"), stringsAsFactors = FALSE) %>%
  as_tibble()

queries <- myco %>%
  mutate(TypeStrain_pident = suppressWarnings(as.numeric(Best_Hit_Pident)),
         TierA_reason = case_when(
           is.na(TypeStrain_pident)            ~ "No type-strain hit",
           TypeStrain_pident < PIDENT_FAMILY   ~ "Identity < 86.5%",
           !grepl("Mycoplasm", Best_Hit_Organism) ~ "Best hit outside Mycoplasmataceae",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(TierA_reason), !is.na(Sequence)) %>%
  dplyr::select(Label, ASV_ID, Sequence, N_Samples, Mean_RelAbund_Pct,
                TypeStrain_hit = Best_Hit_Organism, TypeStrain_pident, TierA_reason)

cat("Queries:", nrow(queries), "\n")
print(queries %>% dplyr::select(Label, TypeStrain_pident, TierA_reason), n = Inf)

fasta <- paste0(">", queries$Label, "\n", queries$Sequence, collapse = "\n")

# =============================================================================
# STEP 2: Submit to NCBI BLAST (URL API), poll, fetch JSON
# =============================================================================
# NCBI guidelines: one submission, poll no more than once per minute.

raw_json <- file.path(blast_dir, "Myco_core_nt_raw.json")

if (!file.exists(raw_json)) {
  put <- POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi", encode = "form", body = list(
    CMD = "Put", PROGRAM = "blastn", MEGABLAST = "on", DATABASE = "core_nt",
    QUERY = fasta, HITLIST_SIZE = "100", EXPECT = "1e-10",
    TOOL = "octocoral_microbiome_revision", EMAIL = Sys.getenv("NCBI_EMAIL")))
  stop_for_status(put)
  txt <- content(put, "text", encoding = "UTF-8")
  rid <- str_match(txt, "RID = (\\S+)")[, 2]
  rtoe <- as.numeric(str_match(txt, "RTOE = (\\d+)")[, 2])
  if (is.na(rid)) stop("No RID returned by NCBI.")
  cat("Submitted. RID =", rid, "| estimated time", rtoe, "s\n")
  writeLines(rid, file.path(blast_dir, "Myco_core_nt_RID.txt"))

  Sys.sleep(max(rtoe, 30))
  repeat {
    st <- content(GET("https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                      query = list(CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid)),
                  "text", encoding = "UTF-8")
    status <- str_match(st, "Status=(\\w+)")[, 2]
    cat(format(Sys.time(), "%H:%M:%S"), "status:", status, "\n")
    if (identical(status, "READY")) break
    if (status %in% c("FAILED", "UNKNOWN")) stop("BLAST search ", status)
    Sys.sleep(60)
  }

  res <- GET("https://blast.ncbi.nlm.nih.gov/Blast.cgi",
             query = list(CMD = "Get", RID = rid, FORMAT_TYPE = "JSON2_S"))
  stop_for_status(res)
  writeLines(content(res, "text", encoding = "UTF-8"), raw_json)
  cat("Saved", raw_json, "\n")
} else {
  cat("Using saved result:", raw_json, "\n")
}

# =============================================================================
# STEP 3: Parse hits
# =============================================================================

bo <- fromJSON(raw_json, simplifyVector = FALSE)$BlastOutput2

hits <- map_dfr(bo, function(x) {
  s <- x$report$results$search
  qlen <- s$query_len
  map_dfr(seq_along(s$hits), function(i) {
    h <- s$hits[[i]]; hsp <- h$hsps[[1]]; desc <- h$description[[1]]
    tibble(Label = s$query_title, Hit_rank = i,
           NCBI_accession = desc$accession, NCBI_title = desc$title,
           NCBI_taxid = as.character(desc$taxid %||% NA), NCBI_sciname = desc$sciname %||% NA_character_,
           NCBI_pident = round(100 * hsp$identity / hsp$align_len, 2),
           NCBI_align_length = hsp$align_len,
           NCBI_query_coverage_pct = round(100 * (hsp$query_to - hsp$query_from + 1) / qlen, 1),
           NCBI_evalue = hsp$evalue, NCBI_bitscore = hsp$bit_score)
  })
})
cat("Hits parsed:", nrow(hits), "for", n_distinct(hits$Label), "ASVs\n")

# =============================================================================
# STEP 4: NCBI lineage for every hit taxid (batched efetch, ranks by name)
# =============================================================================

taxids <- unique(na.omit(hits$NCBI_taxid))
lineage <- map_dfr(split(taxids, ceiling(seq_along(taxids) / 100)), function(ids) {
  Sys.sleep(0.4)
  x <- read_xml(content(GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                            query = list(db = "taxonomy", id = paste(ids, collapse = ","),
                                         retmode = "xml")), "text", encoding = "UTF-8"))
  map_dfr(xml_find_all(x, "/TaxaSet/Taxon"), function(t) {
    ranks <- xml_text(xml_find_all(t, "./LineageEx/Taxon/Rank"))
    names <- xml_text(xml_find_all(t, "./LineageEx/Taxon/ScientificName"))
    own_rank <- xml_text(xml_find_first(t, "./Rank"))
    own_name <- xml_text(xml_find_first(t, "./ScientificName"))
    get <- function(r) if (own_rank == r) own_name else { v <- names[ranks == r]; if (length(v)) v[1] else NA_character_ }
    tibble(NCBI_taxid = xml_text(xml_find_first(t, "./TaxId")),
           NCBI_class = get("class"), NCBI_order = get("order"),
           NCBI_family = get("family"), NCBI_genus = get("genus"))
  })
})

hits <- hits %>% left_join(lineage, by = "NCBI_taxid")

write.csv(hits %>% filter(Hit_rank <= 10),
          file.path(blast_dir, "Myco_core_nt_hits_top10.csv"), row.names = FALSE)

# =============================================================================
# STEP 5: Tier B decision per ASV
# =============================================================================

best_any <- hits %>% group_by(Label) %>% slice_max(NCBI_bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Label, Best_any_hit = NCBI_title, Best_any_pident = NCBI_pident,
            Best_any_coverage = NCBI_query_coverage_pct)

best_classified <- hits %>% filter(!is.na(NCBI_family)) %>%
  group_by(Label) %>% arrange(desc(NCBI_bitscore), desc(NCBI_pident), .by_group = TRUE) %>%
  dplyr::slice(1) %>% ungroup() %>%
  transmute(Label, Best_classified_hit = NCBI_title, Best_classified_family = NCBI_family,
            Best_classified_genus = NCBI_genus, Best_classified_pident = NCBI_pident)

validation <- queries %>%
  dplyr::select(-Sequence) %>%
  left_join(best_any, by = "Label") %>%
  left_join(best_classified, by = "Label") %>%
  mutate(TierB_status = case_when(
    is.na(Best_classified_pident) ~ "Unverified (no classified core_nt hit)",
    Best_classified_family %in% myco_families & Best_classified_pident >= PIDENT_GENUS &
      grepl("Mycoplasma", Best_classified_genus) ~ "Verified (genus) - core_nt",
    Best_classified_family %in% myco_families & Best_classified_pident >= PIDENT_FAMILY ~ "Verified (family) - core_nt",
    !(Best_classified_family %in% myco_families) ~ "Unverified (family discordant)",
    TRUE ~ "Unverified (identity < 86.5%)")) %>%
  arrange(as.numeric(sub("Myco_", "", Label)))

write.csv(validation, file.path(out_dir, "Myco_core_nt_validation.csv"), row.names = FALSE)

cat("\n=== Tier B result ===\n")
print(validation %>% dplyr::select(Label, N_Samples, TypeStrain_pident, Best_any_pident,
                                   Best_classified_family, Best_classified_pident, TierB_status),
      n = Inf, width = 200)
cat("\n"); print(count(validation, TierB_status))
