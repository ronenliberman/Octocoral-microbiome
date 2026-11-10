# =============================================================================
# 01d_ANCOMBC2_typestrain_BLAST.R — Type-strain BLAST of the ASVs with
#                                   significant ANCOM-BC2 results (Appendix 1B)
# =============================================================================
# Appendix 1B previously came from a web BLAST of 13 Mar 2026
# (Output/V8941X58016-Alignment-HitTable.csv) run on an earlier
# differential-abundance result set: only 12 of its 58 ASVs are significant in
# the current analysis (Figure 4 / Table S7). This script BLASTs the current set.
#
#   1. Collects the unique ASVs in ANCOMBC2_PEE2/Figure4_ANCOMBC2_*_significant.csv
#   2. Submits them to NCBI web BLAST (URL API): blastn, megablast, database
#      rRNA_typestrains/16S_ribosomal_RNA (the database used for Appendix 1A).
#      Queries without a megablast hit are re-submitted with blastn (dc-megablast
#      is not available through the URL API).
#   3. Writes out_dir/ANCOMBC2_r3_16S_hits.csv in the 12-column web hit-table
#      format read by 01_BLAST_validation.R (no header).
# Needs internet. Re-running re-submits the search.
# NOTE (7 Oct 2026): NOT used for the final Sheet 1B. The URL API queue stalled and the
# website's multi-query hit-table downloads were incomplete, so the hit table was made with
# 01d_ANCOMBC2_typestrain_BLAST_local.ps1 (local BLAST+, same database). This script still
# writes ANCOMBC2_r3_ASVs.fasta, the query file used by the local runner.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))
suppressPackageStartupMessages(library(httr))
cat("\n========== 01d: ANCOM-BC2 ASVs vs 16S type strains ==========\n\n")

sig_files <- list.files(file.path(out_dir, "ANCOMBC2_PEE2"), "_significant\\.csv$", full.names = TRUE)
sig <- bind_rows(lapply(sig_files, read.csv, stringsAsFactors = FALSE))
ids <- unique(str_extract(sig$Taxa, "[0-9a-f]{32}$"))
rep_seqs <- load_rep_seqs()
stopifnot(all(ids %in% names(rep_seqs)))
cat("Significant results:", nrow(sig), "| unique ASVs:", length(ids), "\n")
q <- rep_seqs[ids]
writeXStringSet(q, file.path(out_dir, "ANCOMBC2_r3_ASVs.fasta"))

BLAST_URL <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
run_blast <- function(seqs, megablast = TRUE) {
  fasta <- paste0(">", names(seqs), "\n", as.character(seqs), collapse = "\n")
  args <- list(CMD = "Put", PROGRAM = "blastn", DATABASE = "rRNA_typestrains/16S_ribosomal_RNA",
               QUERY = fasta, HITLIST_SIZE = "100", TOOL = "octocoral_r3", EMAIL = Sys.getenv("NCBI_EMAIL"))
  if (megablast) args$MEGABLAST <- "on"
  r <- POST(BLAST_URL, body = args, encode = "form")
  txt <- content(r, "text", encoding = "UTF-8")
  rid <- str_match(txt, "RID = (\\S+)")[, 2]
  if (is.na(rid)) stop("no RID returned")
  cat("  RID", rid, if (megablast) "(megablast)" else "(blastn)", "\n")
  repeat {
    Sys.sleep(30)
    s <- content(GET(BLAST_URL, query = list(CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid)), "text", encoding = "UTF-8")
    st <- str_match(s, "Status=(\\w+)")[, 2]
    cat("   status:", st, "\n")
    if (identical(st, "READY")) break
    if (st %in% c("FAILED", "UNKNOWN")) stop("BLAST search ", st)
  }
  tab <- content(GET(BLAST_URL, query = list(CMD = "Get", RID = rid, FORMAT_TYPE = "Tabular",
                                             ALIGNMENTS = "100", DESCRIPTIONS = "100")), "text", encoding = "UTF-8")
  lines <- strsplit(tab, "\n")[[1]]
  lines <- lines[!grepl("^#|^\\s*$|^<", lines)]
  if (!length(lines)) return(NULL)
  h <- read.delim(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  h <- h[, 1:12]
  names(h) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore")
  # subject ids come as gi|...|ref|NR_xxxxxx.x| — keep the accession
  h$sseqid <- ifelse(grepl("ref\\|", h$sseqid), str_match(h$sseqid, "ref\\|([^|]+)")[, 2], h$sseqid)
  h
}

hits <- run_blast(q, megablast = TRUE)
nohit <- setdiff(ids, hits$qseqid)
cat("Queries with megablast hits:", n_distinct(hits$qseqid), "of", length(ids), "\n")
if (length(nohit)) {
  cat("Re-submitting", length(nohit), "queries with blastn\n")
  h2 <- run_blast(q[nohit], megablast = FALSE)
  hits <- bind_rows(hits, h2)
}
cat("Queries with hits:", n_distinct(hits$qseqid), "of", length(ids),
    "| still without a hit:", paste(setdiff(ids, hits$qseqid), collapse = " "), "\n")
stopifnot(all(hits$qseqid %in% ids))
write.table(hits, file.path(blast_dir, "ANCOMBC2_r3_16S_hits.csv"), sep = ",", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
writeLines(c(paste("Searched:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
             "NCBI web BLAST URL API, blastn megablast (blastn for queries without a megablast hit),",
             "database rRNA_typestrains/16S_ribosomal_RNA, hitlist 100",
             paste("Queries:", length(ids), "| with hits:", n_distinct(hits$qseqid))),
           file.path(blast_dir, "ANCOMBC2_r3_16S_hits_SEARCH_INFO.txt"))
cat("Written: ANCOMBC2_r3_16S_hits.csv\n")
