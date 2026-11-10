# =============================================================================
# 07_TableS8_revised.R — Table S8 (pairwise identity among core ASVs) and
#                        Table S9 (validation status of every ASV)
# =============================================================================
# Editor comment addressed: the tree "suggests ASVs from a host species'
# colonies are all nearly identical". R1 Table S8 gave only raw base
# differences. The revised table gives, for each pair of core ASVs:
#   base differences, number of positions compared (both sequences
#   ungapped), and % identity = 100 * (1 - differences / compared).
# Uses the round-3 alignments from 03/04 (same alignments as the trees).
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

cat("\n========== 07: Tables S8 and S9 ==========\n\n")

pairwise_identity <- function(aln_path, keep_labels) {
  aln <- toupper(as.character(read.dna(aln_path, format = "fasta")))
  keep_labels <- intersect(keep_labels, rownames(aln))
  pairs <- t(combn(keep_labels, 2))
  bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    a <- aln[pairs[i, 1], ]; b <- aln[pairs[i, 2], ]
    ok <- a %in% c("A", "C", "G", "T") & b %in% c("A", "C", "G", "T")
    diffs <- sum(a[ok] != b[ok])
    tibble(ASV_1 = pairs[i, 1], ASV_2 = pairs[i, 2],
           Positions_compared = sum(ok),
           Base_differences   = diffs,
           Percent_identity   = round(100 * (1 - diffs / sum(ok)), 2))
  }))
}

to_matrix <- function(pw, labels, value) {
  m <- matrix(NA, length(labels), length(labels), dimnames = list(labels, labels))
  for (i in seq_len(nrow(pw))) {
    m[pw$ASV_1[i], pw$ASV_2[i]] <- pw[[value]][i]
    m[pw$ASV_2[i], pw$ASV_1[i]] <- pw[[value]][i]
  }
  diag(m) <- if (value == "Percent_identity") 100 else 0
  as.data.frame(m) %>% rownames_to_column("ASV")
}

endo_core <- paste0("Endo_", 1:6)
myco_core <- paste0("Myco_", 1:10)

endo_pw <- pairwise_identity(file.path(mac_dir, "Endo_r3_mafft_aligned.fasta"), endo_core) %>%
  mutate(Group = "Endozoicomonadaceae")
myco_pw <- pairwise_identity(file.path(mac_dir, "Myco_r3_mafft_aligned.fasta"), myco_core) %>%
  mutate(Group = "Mycoplasmataceae")

missing_myco <- setdiff(myco_core, c(myco_pw$ASV_1, myco_pw$ASV_2))
if (length(missing_myco) > 0) cat("Not in the Myco alignment:", missing_myco, "\n")

tableS8 <- bind_rows(endo_pw, myco_pw) %>% relocate(Group)

cat("Endo core % identity range:", paste(range(endo_pw$Percent_identity), collapse = "–"), "\n")
cat("Myco core % identity range:", paste(range(myco_pw$Percent_identity), collapse = "–"), "\n")

# Lower triangle = % identity, upper triangle = base differences (one matrix per genus)
combined_matrix <- function(pw, labels) {
  idn <- to_matrix(pw, labels, "Percent_identity")
  bp  <- to_matrix(pw, labels, "Base_differences")
  m <- as.matrix(idn[, -1]); b <- as.matrix(bp[, -1])
  out <- matrix("", nrow(m), ncol(m), dimnames = dimnames(m))
  out[lower.tri(out)] <- sprintf("%.1f", m[lower.tri(m)])
  out[upper.tri(out)] <- as.character(b[upper.tri(b)])
  diag(out) <- "–"
  as.data.frame(out) %>% rownames_to_column("ASV")
}
endo_mat <- combined_matrix(endo_pw, intersect(endo_core, c(endo_pw$ASV_1, endo_pw$ASV_2)))
myco_mat <- combined_matrix(myco_pw, intersect(myco_core, c(myco_pw$ASV_1, myco_pw$ASV_2)))

# ---- Table S9: every Endozoicomonadaceae / Mycoplasmataceae ASV --------------
s9 <- bind_rows(
  read.csv(file.path(out_dir, "Endo_all_ASVs_with_status.csv"), stringsAsFactors = FALSE),
  read.csv(file.path(out_dir, "Myco_all_ASVs_with_status.csv"), stringsAsFactors = FALSE)
) %>%
  transmute(Group = Query_set, Stable_ID, ASV_ID, N_Samples, Total_Reads,
            Mean_RelAbund_Pct, Host_Species, Dominant_Host,
            SILVA_family, SILVA_genus, NCBI_best_hit = NCBI_organism,
            NCBI_pident, Validation_status,
            In_tree = (grepl("^Verified", Validation_status) & N_Samples >= MIN_SAMPLES_TREE) |
              Stable_ID %in% c(endo_core, myco_core)) %>%
  arrange(Group, desc(Total_Reads))

write.csv(tableS8, file.path(out_dir, "TableS8_core_ASV_pairwise_identity.csv"), row.names = FALSE)
write.csv(s9, file.path(out_dir, "TableS9_ASV_validation_all.csv"), row.names = FALSE)
write_xlsx(list(S8_pairs = tableS8,
                S8_Endo_matrix = endo_mat,
                S8_Myco_matrix = myco_mat,
                S9_all_ASVs = s9),
           file.path(out_dir, "TablesS8_S9_revised.xlsx"))

cat("Written: TableS8_core_ASV_pairwise_identity.csv, TableS9_ASV_validation_all.csv, TablesS8_S9_revised.xlsx\n")
cat("S8 matrix layout: lower triangle = % identity, upper triangle = base differences.\n")
