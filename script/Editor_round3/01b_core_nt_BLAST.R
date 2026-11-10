# =============================================================================
# 01b_core_nt_BLAST.R — Remote BLAST for the ASVs that script 01 could not verify
# =============================================================================
# Run AFTER 01 (first pass). Uses 01's own outputs:
#   out_dir/Unverified_ASVs_for_nt_BLAST.fasta   (Endo + Myco ASVs not verified)
#   out_dir/ASV_validation_status.csv
#
# Two searches, submitted to the NCBI BLAST URL API (blastn, megablast):
#   1. Tier B: every ASV in the FASTA vs core_nt (includes uncultured clones).
#   2. Tier A top-up: ASVs 01 marked "Not yet BLASTed" (present in the PEE2
#      table but not in the R1 BLAST files) vs the 16S rRNA type-strain
#      database, so they get the same first check as all other ASVs.
#
# Writes files that 01 reads automatically on its next run:
#   Unverified_ASVs_nt_hits.csv            core_nt hits WITH family-level taxonomy
#                                          (top 5 per ASV; 12-column hit table)
#   TierA_topup_16S_hits_Endozoicomonadaceae.csv / _Mycoplasmataceae.csv
# plus, for reading:
#   core_nt_raw.json, core_nt_hits_top10.csv (all hits incl. "uncultured bacterium")
#   core_nt_TierB_summary.csv              one row per ASV, Tier B decision
#
# Raw JSON is cached: re-running does not resubmit unless the file is deleted.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(xml2)
  library(Biostrings)
})

rev_dir   <- here::here()
out_dir   <- file.path(rev_dir, "output", "Editor_round3")
blast_dir <- file.path(rev_dir, "blast")   # NCBI results are saved here and reused

PIDENT_FAMILY <- 86.5
PIDENT_GENUS  <- 94.5
family_equivalents <- list(
  Endozoicomonadaceae = c("Endozoicomonadaceae"),
  Mycoplasmataceae    = c("Mycoplasmataceae", "Metamycoplasmataceae", "Mycoplasmoidaceae")
)

cat("\n========== 01b: remote BLAST of unverified ASVs ==========\n\n")

# =============================================================================
# STEP 1: Queries from script 01
# =============================================================================

seqs <- readDNAStringSet(file.path(out_dir, "Unverified_ASVs_for_nt_BLAST.fasta"))
validation <- read.csv(file.path(out_dir, "ASV_validation_status.csv"), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  filter(Query_set %in% names(family_equivalents)) %>%
  distinct(ASV_ID, .keep_all = TRUE)

queries <- tibble(ASV_ID = names(seqs), Sequence = as.character(seqs)) %>%
  left_join(validation %>% dplyr::select(ASV_ID, Query_set, SILVA_family, SILVA_genus,
                                         TierA_status = Validation_status,
                                         TierA_hit = NCBI_organism, TierA_pident = NCBI_pident),
            by = "ASV_ID")
cat("Queries:", nrow(queries), "\n")
print(count(queries, Query_set, TierA_status))

# =============================================================================
# STEP 2: BLAST URL API helper (one submission per database, cached)
# =============================================================================

run_remote_blast <- function(q, database, cache_json) {
  if (file.exists(cache_json)) {
    cat("Using cached", basename(cache_json), "\n")
    return(fromJSON(cache_json, simplifyVector = FALSE)$BlastOutput2)
  }
  fasta <- paste0(">", q$ASV_ID, "\n", q$Sequence, collapse = "\n")
  put <- POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi", encode = "form", body = list(
    CMD = "Put", PROGRAM = "blastn", MEGABLAST = "on", DATABASE = database,
    QUERY = fasta, HITLIST_SIZE = "100", EXPECT = "1e-10",
    TOOL = "octocoral_microbiome_revision", EMAIL = Sys.getenv("NCBI_EMAIL")))
  stop_for_status(put)
  txt <- content(put, "text", encoding = "UTF-8")
  rid <- str_match(txt, "RID = (\\S+)")[, 2]
  rtoe <- as.numeric(str_match(txt, "RTOE = (\\d+)")[, 2])
  if (is.na(rid)) stop("No RID returned by NCBI.")
  cat(database, ": RID", rid, "| estimated", rtoe, "s\n")
  write(paste(database, rid), file.path(blast_dir, "BLAST_RIDs.txt"), append = TRUE)

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
  writeLines(content(res, "text", encoding = "UTF-8"), cache_json)
  fromJSON(cache_json, simplifyVector = FALSE)$BlastOutput2
}

parse_hits <- function(bo) {
  map_dfr(bo, function(x) {
    s <- x$report$results$search
    qlen <- s$query_len
    map_dfr(seq_along(s$hits), function(i) {
      h <- s$hits[[i]]; hsp <- h$hsps[[1]]; desc <- h$description[[1]]
      tibble(ASV_ID = s$query_title, Hit_rank = i,
             NCBI_accession = desc$accession, NCBI_title = desc$title,
             NCBI_taxid = as.character(desc$taxid %||% NA),
             pident = round(100 * hsp$identity / hsp$align_len, 3),
             length = hsp$align_len, mismatch = hsp$align_len - hsp$identity - (hsp$gaps %||% 0),
             gapopen = hsp$gaps %||% 0, qstart = hsp$query_from, qend = hsp$query_to,
             sstart = hsp$hit_from, send = hsp$hit_to, evalue = hsp$evalue,
             bitscore = hsp$bit_score,
             coverage_pct = round(100 * (hsp$query_to - hsp$query_from + 1) / qlen, 1))
    })
  })
}

# Lineage for taxids, batched (100 per request), ranks looked up by name
fetch_lineages <- function(taxids) {
  taxids <- unique(na.omit(taxids))
  map_dfr(split(taxids, ceiling(seq_along(taxids) / 100)), function(ids) {
    Sys.sleep(0.4)
    x <- read_xml(content(GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                              query = list(db = "taxonomy", id = paste(ids, collapse = ","),
                                           retmode = "xml")), "text", encoding = "UTF-8"))
    map_dfr(xml_find_all(x, "/TaxaSet/Taxon"), function(t) {
      ranks <- xml_text(xml_find_all(t, "./LineageEx/Taxon/Rank"))
      nms   <- xml_text(xml_find_all(t, "./LineageEx/Taxon/ScientificName"))
      own_rank <- xml_text(xml_find_first(t, "./Rank"))
      own_name <- xml_text(xml_find_first(t, "./ScientificName"))
      get <- function(r) if (own_rank == r) own_name else { v <- nms[ranks == r]; if (length(v)) v[1] else NA_character_ }
      tibble(NCBI_taxid = xml_text(xml_find_first(t, "./TaxId")), NCBI_sciname = own_name,
             NCBI_family = get("family"), NCBI_genus = get("genus"))
    })
  })
}

hit_table_cols <- c("ASV_ID", "NCBI_accession", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# =============================================================================
# STEP 3: Tier B — core_nt
# =============================================================================

nt_hits <- parse_hits(run_remote_blast(queries, "core_nt", file.path(blast_dir, "core_nt_raw.json")))
nt_hits <- nt_hits %>% left_join(fetch_lineages(nt_hits$NCBI_taxid), by = "NCBI_taxid")
cat("core_nt hits:", nrow(nt_hits), "for", n_distinct(nt_hits$ASV_ID), "ASVs\n")

write.csv(nt_hits %>% filter(Hit_rank <= 10) %>% left_join(queries %>% dplyr::select(ASV_ID, Query_set), by = "ASV_ID"),
          file.path(blast_dir, "core_nt_hits_top10.csv"), row.names = FALSE)

# Hit table for 01: only hits that carry a family, top 5 per ASV
nt_classified <- nt_hits %>%
  filter(!is.na(NCBI_family)) %>%
  group_by(ASV_ID) %>% arrange(desc(bitscore), desc(pident), .by_group = TRUE) %>%
  slice_head(n = 5) %>% ungroup()
write.table(nt_classified[, hit_table_cols], file.path(blast_dir, "Unverified_ASVs_nt_hits.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# =============================================================================
# STEP 4: Tier A top-up for ASVs never BLASTed against the type-strain database
# =============================================================================

topup <- queries %>% filter(TierA_status == "Not yet BLASTed")
if (nrow(topup) > 0) {
  cat("\nTier A top-up (16S type strains):", nrow(topup), "ASVs\n")
  ts_hits <- parse_hits(run_remote_blast(topup, "rRNA_typestrains/16S_ribosomal_RNA",
                                         file.path(blast_dir, "typestrain_topup_raw.json")))
  ts_hits <- ts_hits %>% left_join(topup %>% dplyr::select(ASV_ID, Query_set), by = "ASV_ID")
  for (set in unique(ts_hits$Query_set)) {
    write.table(ts_hits %>% filter(Query_set == set, Hit_rank <= 10) %>% dplyr::select(all_of(hit_table_cols)),
                file.path(blast_dir, paste0("TierA_topup_16S_hits_", set, ".csv")),
                sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  cat("Top-up hits written for:", paste(unique(ts_hits$Query_set), collapse = ", "), "\n")
}

# =============================================================================
# STEP 5: Tier B decision (for reading; 01 applies the same rule)
# =============================================================================

best_any <- nt_hits %>% group_by(ASV_ID) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>% transmute(ASV_ID, Best_any_hit = NCBI_title, Best_any_pident = pident,
                          Best_any_coverage = coverage_pct)
best_cls <- nt_classified %>% group_by(ASV_ID) %>% dplyr::slice(1) %>% ungroup() %>%
  transmute(ASV_ID, Best_classified_hit = NCBI_title, Best_classified_family = NCBI_family,
            Best_classified_genus = NCBI_genus, Best_classified_pident = pident,
            Best_classified_accession = NCBI_accession)

summary_tbl <- queries %>%
  dplyr::select(-Sequence) %>%
  left_join(best_any, by = "ASV_ID") %>%
  left_join(best_cls, by = "ASV_ID") %>%
  mutate(
    fam_ok = map2_lgl(Query_set, Best_classified_family, ~ !is.na(.y) && .y %in% family_equivalents[[.x]]),
    TierB_status = case_when(
      is.na(Best_classified_pident) ~ "Unverified (no classified core_nt hit)",
      fam_ok & Best_classified_pident >= PIDENT_GENUS & Best_classified_genus == SILVA_genus ~ "Verified (genus) - core_nt",
      fam_ok & Best_classified_pident >= PIDENT_FAMILY ~ "Verified (family) - core_nt",
      !fam_ok ~ "Unverified (family discordant)",
      TRUE ~ "Unverified (identity < 86.5%)")) %>%
  dplyr::select(-fam_ok) %>%
  arrange(Query_set, TierB_status)

write.csv(summary_tbl, file.path(out_dir, "core_nt_TierB_summary.csv"), row.names = FALSE)

cat("\n=== Tier B summary ===\n")
print(count(summary_tbl, Query_set, TierB_status), n = Inf)
cat("\nNext: re-run 01_BLAST_validation.R (uses the cache; folds in Tier A top-up and Tier B), then 02.\n")
