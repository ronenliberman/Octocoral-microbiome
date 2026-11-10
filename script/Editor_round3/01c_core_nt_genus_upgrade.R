# =============================================================================
# 01c_core_nt_genus_upgrade.R — Can family-level ASVs be named at genus level?
# =============================================================================
# ASVs that 01 verified only at FAMILY level (same family as SILVA, 86.5–94.5%
# identity to the closest type strain, or a different genus) are BLASTed
# against core_nt. On its next run, 01 upgrades an ASV to "Verified (genus)"
# when its best core_nt hit WITH genus-level taxonomy is:
#   same family as SILVA, same genus as SILVA, and >= 94.5% identity.
# Nothing is ever downgraded by this step.
#
# Run order: 01 -> 01b -> 01 -> 01c -> 01 -> 02
#
# Input : out_dir/ASV_validation_status.csv (from 01)
# Output: out_dir/Family_level_ASVs_for_nt_BLAST.fasta
#         out_dir/genus_check_core_nt_raw.json      (cached; delete to resubmit)
#         out_dir/Family_level_ASVs_nt_hits.csv     (12-col hit table, read by 01)
#         out_dir/genus_check_core_nt_hits_top10.csv (all hits, for reading)
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
raw_dir   <- file.path(rev_dir, "raw files")

cat("\n========== 01c: core_nt genus check for family-level ASVs ==========\n\n")

validation <- read.csv(file.path(out_dir, "ASV_validation_status.csv"), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  filter(Query_set %in% c("Endozoicomonadaceae", "Mycoplasmataceae"),
         Validation_status == "Verified (family)") %>%
  distinct(ASV_ID, .keep_all = TRUE)

rep_seqs <- readDNAStringSet(file.path(raw_dir, "dna-sequences-PEE2.fasta"))
queries <- validation %>%
  transmute(ASV_ID, Query_set, SILVA_genus, TypeStrain_hit = NCBI_organism,
            TypeStrain_pident = NCBI_pident,
            Sequence = as.character(rep_seqs[ASV_ID]))
cat("Family-level ASVs to test:", nrow(queries), "\n")
print(count(queries, Query_set))

qs <- DNAStringSet(queries$Sequence); names(qs) <- queries$ASV_ID
writeXStringSet(qs, file.path(out_dir, "Family_level_ASVs_for_nt_BLAST.fasta"))

# ---- BLAST results: website download (preferred) or remote URL API ----------
# WEBSITE MODE: BLAST Family_level_ASVs_for_nt_BLAST.fasta on the NCBI website
# (database core_nt, megablast, 100 target sequences), download
# "Hit Table (CSV)" and save it as blast/genus_check_core_nt_WEB_HitTable.csv.
# If that file exists, the remote search below is skipped.
# Several downloads are allowed (the website can drop queries from large jobs):
# every file named genus_check_core_nt_WEB_HitTable*.csv is read and combined.
web_files  <- list.files(blast_dir, pattern = "^genus_check_core_nt_WEB_HitTable.*\\.csv$", full.names = TRUE)
cache_json <- file.path(blast_dir, "genus_check_core_nt_raw.json")

if (length(web_files) > 0) {
  cat("Website mode: reading", paste(basename(web_files), collapse = ", "), "\n")
  qlen <- setNames(nchar(queries$Sequence), queries$ASV_ID)
  hits <- map_dfr(web_files, ~ read.csv(.x, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("qseqid", "NCBI_accession", "pident", "length", "mismatch",
                                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))) %>%
    as_tibble() %>%
    mutate(ASV_ID = str_extract(qseqid, "[0-9a-f]{32}")) %>%
    filter(!is.na(ASV_ID)) %>%
    group_by(ASV_ID) %>%
    arrange(desc(bitscore), desc(pident), .by_group = TRUE) %>%
    distinct(NCBI_accession, .keep_all = TRUE) %>%
    mutate(Hit_rank = row_number()) %>%
    ungroup() %>%
    mutate(coverage_pct = round(100 * (qend - qstart + 1) / qlen[ASV_ID], 1),
           NCBI_title = NA_character_)
  cat("Web hits:", nrow(hits), "for", n_distinct(hits$ASV_ID), "ASVs\n")
  missing_q <- setdiff(queries$ASV_ID, hits$ASV_ID)
  if (length(missing_q) > 0) {
    cat("WARNING:", length(missing_q), "of", nrow(queries),
        "queries have NO hits in the web files — BLAST them again (genus_check_batch*.fasta):\n")
    print(queries %>% filter(ASV_ID %in% missing_q) %>% dplyr::select(Query_set, ASV_ID, TypeStrain_hit))
  }

  # accession -> taxid (+ title), batched esummary (200 per request)
  accs <- unique(hits$NCBI_accession)
  acc_tax <- map_dfr(split(accs, ceiling(seq_along(accs) / 200)), function(ids) {
    Sys.sleep(0.4)
    js <- fromJSON(content(POST("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                                encode = "form",
                                body = list(db = "nuccore", id = paste(ids, collapse = ","), retmode = "json")),
                           "text", encoding = "UTF-8"), simplifyVector = FALSE)$result
    uids <- setdiff(names(js), "uids")
    map_dfr(uids, function(u) tibble(NCBI_accession = js[[u]]$accessionversion,
                                     NCBI_taxid = as.character(js[[u]]$taxid),
                                     title = js[[u]]$title))
  })
  hits <- hits %>%
    left_join(acc_tax, by = "NCBI_accession") %>%
    mutate(NCBI_title = title) %>%
    dplyr::select(-title, -qseqid)
  cat("Accessions with taxid:", sum(!is.na(acc_tax$NCBI_taxid)), "of", length(accs), "\n")

} else {
if (!file.exists(cache_json)) {
  fasta <- paste0(">", queries$ASV_ID, "\n", queries$Sequence, collapse = "\n")
  put <- POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi", encode = "form", body = list(
    CMD = "Put", PROGRAM = "blastn", MEGABLAST = "on", DATABASE = "core_nt",
    QUERY = fasta, HITLIST_SIZE = "100", EXPECT = "1e-10",
    TOOL = "octocoral_microbiome_revision", EMAIL = Sys.getenv("NCBI_EMAIL")))
  stop_for_status(put)
  txt <- content(put, "text", encoding = "UTF-8")
  rid <- str_match(txt, "RID = (\\S+)")[, 2]
  if (is.na(rid)) stop("No RID returned by NCBI.")
  cat("core_nt: RID", rid, "\n")
  write(paste("core_nt genus check", rid), file.path(blast_dir, "BLAST_RIDs.txt"), append = TRUE)
  Sys.sleep(30)
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
} else {
  cat("Using cached", basename(cache_json), "\n")
}

bo <- fromJSON(cache_json, simplifyVector = FALSE)$BlastOutput2
hits <- map_dfr(bo, function(x) {
  s <- x$report$results$search
  map_dfr(seq_along(s$hits), function(i) {
    h <- s$hits[[i]]; hsp <- h$hsps[[1]]; desc <- h$description[[1]]
    tibble(ASV_ID = s$query_title, Hit_rank = i, NCBI_accession = desc$accession,
           NCBI_title = desc$title, NCBI_taxid = as.character(desc$taxid %||% NA),
           pident = round(100 * hsp$identity / hsp$align_len, 3), length = hsp$align_len,
           mismatch = hsp$align_len - hsp$identity - (hsp$gaps %||% 0), gapopen = hsp$gaps %||% 0,
           qstart = hsp$query_from, qend = hsp$query_to, sstart = hsp$hit_from, send = hsp$hit_to,
           evalue = hsp$evalue, bitscore = hsp$bit_score,
           coverage_pct = round(100 * (hsp$query_to - hsp$query_from + 1) / s$query_len, 1))
  })
})
}  # end remote mode

# Lineages (batched), ranks by name
taxids <- unique(na.omit(hits$NCBI_taxid))
lin <- map_dfr(split(taxids, ceiling(seq_along(taxids) / 100)), function(ids) {
  Sys.sleep(0.4)
  x <- read_xml(content(GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                            query = list(db = "taxonomy", id = paste(ids, collapse = ","), retmode = "xml")),
                        "text", encoding = "UTF-8"))
  map_dfr(xml_find_all(x, "/TaxaSet/Taxon"), function(t) {
    ranks <- xml_text(xml_find_all(t, "./LineageEx/Taxon/Rank"))
    nms   <- xml_text(xml_find_all(t, "./LineageEx/Taxon/ScientificName"))
    own_rank <- xml_text(xml_find_first(t, "./Rank")); own_name <- xml_text(xml_find_first(t, "./ScientificName"))
    get <- function(r) if (own_rank == r) own_name else { v <- nms[ranks == r]; if (length(v)) v[1] else NA_character_ }
    tibble(NCBI_taxid = xml_text(xml_find_first(t, "./TaxId")),
           NCBI_family = get("family"), NCBI_genus = get("genus"))
  })
})
hits <- hits %>% left_join(lin, by = "NCBI_taxid")

write.csv(hits %>% filter(Hit_rank <= 10), file.path(blast_dir, "genus_check_core_nt_hits_top10.csv"), row.names = FALSE)

# Hit table for 01: hits with GENUS-level taxonomy, top 5 per ASV
cols <- c("ASV_ID", "NCBI_accession", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore")
genus_hits <- hits %>% filter(!is.na(NCBI_genus)) %>%
  group_by(ASV_ID) %>% arrange(desc(bitscore), desc(pident), .by_group = TRUE) %>%
  slice_head(n = 5) %>% ungroup()
write.table(genus_hits[, cols], file.path(blast_dir, "Family_level_ASVs_nt_hits.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Preview of the decision (01 makes the real one)
preview <- genus_hits %>% group_by(ASV_ID) %>% dplyr::slice(1) %>% ungroup() %>%
  dplyr::select(ASV_ID, nt_genus = NCBI_genus, nt_family = NCBI_family, nt_pident = pident, nt_hit = NCBI_title) %>%
  right_join(queries %>% dplyr::select(-Sequence), by = "ASV_ID") %>%
  mutate(Upgrade = !is.na(nt_pident) & nt_pident >= 94.5 & nt_genus == SILVA_genus)
cat("\nWould upgrade to genus:", sum(preview$Upgrade), "of", nrow(preview), "\n")
print(count(preview, Query_set, Upgrade))
cat("\nNext: re-run 01_BLAST_validation.R, then 02.\n")
