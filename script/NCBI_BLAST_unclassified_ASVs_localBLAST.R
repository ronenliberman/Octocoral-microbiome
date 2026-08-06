################################################################################
# NCBI BLAST Validation of Poorly Classified ASVs - LOCAL BLAST+ VERSION
# Author: Ronen Liberman
# Date: March 2026
# Description: Identifies ASVs not classified above Class level by SILVA 138.2,
#   BLASTs them against NCBI nt using local BLAST+, and produces a CSV with
#   top 3 BLAST hits + NCBI taxonomy alongside SILVA classifications.
#
# This version uses NCBI BLAST+ web service via the -remote flag (no local
# database download needed) OR local nt database if available.
#
# PREREQUISITES:
#   - BLAST+ installed: conda install -c bioconda blast
#   - R packages: tidyverse, microeco, file2meco, qiime2R, Biostrings, rentrez
#
# TWO MODES:
#   Mode 1 (remote, no setup): blastn -remote (slower, ~1-2 min/ASV)
#   Mode 2 (local, fast): requires nt database (~100GB)
################################################################################

# =============================================================================
# SETUP
# =============================================================================

library(tidyverse)
library(microeco)
library(file2meco)
library(qiime2R)
library(Biostrings)
#install.packages("rentrez")
library(rentrez)     # for fetching taxonomy from TaxIDs
library(XML)

setwd("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/")
output_dir <- file.path("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# CONFIGURATION: Set to TRUE if you have a local nt database, FALSE for remote
USE_LOCAL_DB <- FALSE
# If local, set the path to the nt database
LOCAL_DB_PATH <- "/path/to/blastdb/nt"  # change this
# Number of threads for local BLAST
BLAST_THREADS <- 4

# =============================================================================
# SECTION 1: LOAD DATA AND IDENTIFY POORLY CLASSIFIED ASVs
# =============================================================================

cat("=== Loading data ===\n")

data <- qiime2meco(
  "raw_data/table-NoEMC-PEE2.qza",
  taxonomy_table = "raw_data/taxonomy-PEE2.qza",
  sample_table   = "raw_data/metadata_all_edited.csv",
  auto_tidy      = TRUE
)

rep_seqs <- readDNAStringSet("raw_data/dna-sequences-PEE2.fasta")

d <- clone(data)
d$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

# Filter to focal species and sites
d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE", "SW"))
d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

# Remove unreliable samples identified during QC
samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024")

cat("=== REMOVING UNRELIABLE SAMPLES ===\n")
cat("Samples before filtering:", nrow(d$sample_table), "\n")

d$sample_table <- d$sample_table %>%
  filter(!rownames(.) %in% samples_to_remove)

cat("Samples after filtering:", nrow(d$sample_table), "\n\n")

# Rename species to full scientific names
species.names <- list(
  "SW" = "Seawater",
  "BAST" = "B. asbestinum",
  "ECAR" = "E. caribaeorum",
  "MMUR" = "M. muricata",
  "EFLE" = "E. flexuosa"
)

d$sample_table <- d$sample_table %>%
  mutate(SPP = recode(SPP, !!!species.names))

# Tidy dataset to sync all tables
d$tidy_dataset()

cat("Final sample count:", nrow(d$sample_table), "\n")
cat("Species included:", paste(unique(d$sample_table$SPP), collapse = ", "), "\n")

cat("Total ASVs:", nrow(d$tax_table), "\n")
cat("Total rep seqs:", length(rep_seqs), "\n")

# Build taxonomy dataframe
tax_df <- d$tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID")

# Identify poorly classified ASVs: no informative Class-level assignment
# Check for empty, prefix-only, or "unclassified" patterns at Class level
poorly_classified <- tax_df %>%
  filter(
    is.na(Class) | Class == "" | Class == "c__" |
    grepl("^c__$", Class) |
    grepl("^(c__)?unclassified|^(c__)?uncultured|^(c__)?unknown|Incertae",
          Class, ignore.case = TRUE)
  )

cat("Poorly classified ASVs (no Class-level ID):", nrow(poorly_classified), "\n")

# Classification depth summary
cat("\nFull classification depth summary:\n")
depth_summary <- tax_df %>%
  mutate(
    deepest = case_when(
      !(is.na(Genus) | Genus == "" | Genus == "g__" |
        grepl("unclassified|uncultured|unknown", Genus, ignore.case = TRUE)) ~ "Genus",
      !(is.na(Family) | Family == "" | Family == "f__" |
        grepl("unclassified|uncultured|unknown", Family, ignore.case = TRUE)) ~ "Family",
      !(is.na(Order) | Order == "" | Order == "o__" |
        grepl("unclassified|uncultured|unknown", Order, ignore.case = TRUE)) ~ "Order",
      !(is.na(Class) | Class == "" | Class == "c__" |
        grepl("unclassified|uncultured|unknown", Class, ignore.case = TRUE)) ~ "Class",
      !(is.na(Phylum) | Phylum == "" | Phylum == "p__" |
        grepl("unclassified|uncultured|unknown", Phylum, ignore.case = TRUE)) ~ "Phylum",
      TRUE ~ "Kingdom_only"
    )
  )
print(table(depth_summary$deepest))

# Get abundance info for poorly classified ASVs
asvs_to_blast <- poorly_classified$ASV_ID
asv_abundance <- data.frame(
  ASV_ID      = asvs_to_blast,
  total_reads = rowSums(d$otu_table[asvs_to_blast, , drop = FALSE]),
  n_samples   = rowSums(d$otu_table[asvs_to_blast, , drop = FALSE] > 0)
) %>%
  arrange(desc(total_reads))

cat("\nAbundance of poorly classified ASVs:\n")
cat("  Total reads:", format(sum(asv_abundance$total_reads), big.mark = ","), "\n")
cat("  Number of ASVs:", nrow(asv_abundance), "\n")
cat("  Top 5:\n")
print(head(asv_abundance, 5))

# =============================================================================
# SECTION 2: EXTRACT SEQUENCES AND RUN BLAST
# =============================================================================

cat("\n=== Extracting sequences and running BLAST ===\n")

# Extract and save FASTA
blast_seqs <- rep_seqs[names(rep_seqs) %in% asvs_to_blast]
fasta_path <- file.path(output_dir, "poorly_classified_ASVs.fasta")
writeXStringSet(blast_seqs, fasta_path)
cat("FASTA written:", length(blast_seqs), "sequences to", fasta_path, "\n")

# Define BLAST output path
blast_out_path <- file.path(output_dir, "blast_results.tsv")

# Read the NCBI web BLAST CSV download
# The hit table CSV has these columns:
# query acc., subject acc., % identity, alignment length, mismatches,
# gap opens, q. start, q. end, s. start, s. end, evalue, bit score

blast_results <- read.csv(
  file.path(output_dir, "blast_results.csv"),
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

# Rename columns to match our script expectations
colnames(blast_results) <- c("ASV_ID", "sseqid", "pident", "length",
                             "mismatch", "gapopen", "qstart", "qend",
                             "sstart", "send", "evalue", "bitscore")

# Add empty columns for fields not in the web download
blast_results$staxids   <- NA
blast_results$sscinames <- NA
blast_results$stitle    <- NA
head(blast_results)

cat("Total BLAST hits:", nrow(blast_results), "\n")
cat("Unique ASVs with hits:", length(unique(blast_results$ASV_ID)), "\n")

# =============================================================================
# SECTION 3: PARSE BLAST RESULTS
# =============================================================================

cat("\n=== Parsing BLAST results ===\n")

blast_results <- read.delim(
  blast_out_path,
  header = FALSE,
  col.names = c("ASV_ID", "sseqid", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore", "staxids", "sscinames", "stitle"),
  stringsAsFactors = FALSE
)

cat("Total BLAST hits:", nrow(blast_results), "\n")
cat("Unique ASVs with hits:", length(unique(blast_results$ASV_ID)), "\n")

# Rank hits by pident within each ASV, keep top 3
blast_top3 <- blast_results %>%
  group_by(ASV_ID) %>%
  arrange(desc(pident), evalue) %>%
  mutate(hit_rank = row_number()) %>%
  filter(hit_rank <= 3) %>%
  ungroup()

cat("Top 3 hits per ASV:", nrow(blast_top3), "\n")

# =============================================================================
# SECTION 4: FETCH FULL NCBI TAXONOMY FOR BEST HITS
# =============================================================================

cat("\n=== SECTION 4: Fetching NCBI taxonomy for best hits ===\n")

# Get best hit per ASV (highest pident)
best_hits <- blast_top3 %>%
  group_by(ASV_ID) %>%
  slice_max(pident, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Best hits to resolve taxonomy for:", nrow(best_hits), "\n")
cat("Pident range of best hits:", min(best_hits$pident), "-", max(best_hits$pident), "\n")

# Function to fetch full NCBI taxonomy from an accession number
# Works in two steps: accession -> TaxID -> full lineage
fetch_taxonomy_from_accession <- function(accession) {
  tryCatch({
    # Step 1: Get TaxID from the nucleotide record
    nuc_summary <- entrez_summary(db = "nucleotide", id = accession)
    taxid <- nuc_summary$taxid
    
    if (is.null(taxid) || is.na(taxid)) {
      return(data.frame(
        accession = accession, ncbi_taxid = NA,
        ncbi_Kingdom = NA, ncbi_Phylum = NA, ncbi_Class = NA,
        ncbi_Order = NA, ncbi_Family = NA, ncbi_Genus = NA,
        ncbi_Species = NA, stringsAsFactors = FALSE
      ))
    }
    
    # Step 2: Get full taxonomy lineage from TaxID
    tax_record <- entrez_fetch(db = "taxonomy", id = taxid, rettype = "xml")
    tax_xml <- xmlParse(tax_record, asText = TRUE)
    
    # Extract all ranks from the lineage
    lineage_nodes <- getNodeSet(tax_xml, "//LineageEx/Taxon")
    lineage <- list()
    for (node in lineage_nodes) {
      rank <- xmlValue(getNodeSet(node, ".//Rank")[[1]])
      name <- xmlValue(getNodeSet(node, ".//ScientificName")[[1]])
      lineage[[rank]] <- name
    }
    
    # Get the organism-level scientific name
    sci_name <- xmlValue(getNodeSet(tax_xml, "//ScientificName")[[1]])
    free(tax_xml)
    
    data.frame(
      accession    = accession,
      ncbi_taxid   = taxid,
      ncbi_Kingdom = lineage[["superkingdom"]] %||% NA_character_,
      ncbi_Phylum  = lineage[["phylum"]] %||% NA_character_,
      ncbi_Class   = lineage[["class"]] %||% NA_character_,
      ncbi_Order   = lineage[["order"]] %||% NA_character_,
      ncbi_Family  = lineage[["family"]] %||% NA_character_,
      ncbi_Genus   = lineage[["genus"]] %||% NA_character_,
      ncbi_Species = sci_name,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(" FAILED:", conditionMessage(e), "\n")
    data.frame(
      accession = accession, ncbi_taxid = NA,
      ncbi_Kingdom = NA, ncbi_Phylum = NA, ncbi_Class = NA,
      ncbi_Order = NA, ncbi_Family = NA, ncbi_Genus = NA,
      ncbi_Species = NA, stringsAsFactors = FALSE
    )
  })
}

# Fetch taxonomy for each unique accession from best hits
unique_accessions <- unique(best_hits$sseqid)
cat("Unique accessions to look up:", length(unique_accessions), "\n\n")

acc_results <- list()
for (i in seq_along(unique_accessions)) {
  acc <- unique_accessions[i]
  cat(sprintf("  [%d/%d] %s ...", i, length(unique_accessions), acc))
  acc_results[[acc]] <- fetch_taxonomy_from_accession(acc)
  cat(" done\n")
  Sys.sleep(0.5)  # respect NCBI rate limits (3 requests/sec without API key)
}

ncbi_tax_lookup <- bind_rows(acc_results)

cat("\nTaxonomy fetch complete.\n")
cat("Successful lookups:", sum(!is.na(ncbi_tax_lookup$ncbi_taxid)), "/",
    nrow(ncbi_tax_lookup), "\n")

# Join taxonomy back to best hits
best_hits_with_tax <- best_hits %>%
  left_join(ncbi_tax_lookup, by = c("sseqid" = "accession"))

cat("ASVs with NCBI Class resolved:", 
    sum(!is.na(best_hits_with_tax$ncbi_Class) & best_hits_with_tax$ncbi_Class != ""), "\n")
cat("ASVs with NCBI Genus resolved:", 
    sum(!is.na(best_hits_with_tax$ncbi_Genus) & best_hits_with_tax$ncbi_Genus != ""), "\n")
# =============================================================================
# SECTION 5: BUILD FINAL MERGED CSV
# =============================================================================

cat("\n=== Building final output ===\n")

# SILVA taxonomy for poorly classified ASVs
silva_tax <- tax_df %>%
  filter(ASV_ID %in% asvs_to_blast) %>%
  select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  rename_with(~ paste0("silva_", .), -ASV_ID)

# Top 3 hits in wide format
top3_wide <- blast_top3 %>%
  select(ASV_ID, hit_rank, sseqid, sscinames, pident, evalue, stitle) %>%
  pivot_wider(
    id_cols     = ASV_ID,
    names_from  = hit_rank,
    values_from = c(sseqid, sscinames, pident, evalue, stitle),
    names_glue  = "NCBI_hit{hit_rank}_{.value}"
  )

ncbi_best <- best_hits_with_tax %>%
  dplyr::select(ASV_ID, pident, ncbi_taxid, ncbi_Kingdom, ncbi_Phylum, ncbi_Class,
                ncbi_Order, ncbi_Family, ncbi_Genus, ncbi_Species) %>%
  dplyr::rename(best_hit_pident = pident)

# Merge all
final <- silva_tax %>%
  left_join(asv_abundance, by = "ASV_ID") %>%
  left_join(ncbi_best, by = "ASV_ID") %>%
  left_join(top3_wide, by = "ASV_ID") %>%
  arrange(desc(total_reads))

# Save
final_path <- file.path(output_dir, "NCBI_BLAST_poorly_classified_ASVs_summary.csv")
write.csv(final, final_path, row.names = FALSE)

cat("\n========================================\n")
cat("FINAL OUTPUT:", final_path, "\n")
cat("========================================\n")
cat("Total poorly classified ASVs:", nrow(silva_tax), "\n")
cat("ASVs with BLAST hits:", sum(!is.na(final$NCBI_hit1_sseqid)), "\n")
cat("ASVs with NCBI Class resolved:", sum(!is.na(final$ncbi_Class) & final$ncbi_Class != ""), "\n")
cat("ASVs with NCBI Genus resolved:", sum(!is.na(final$ncbi_Genus) & final$ncbi_Genus != ""), "\n")

# Quick comparison
cat("\nImproved classification via NCBI:\n")
improved <- final %>%
  filter(!is.na(ncbi_Class) & ncbi_Class != "") %>%
  select(ASV_ID, total_reads, n_samples, silva_Phylum, ncbi_Class, ncbi_Order,
         ncbi_Family, ncbi_Genus, ncbi_Species, best_hit_pident)
print(head(improved, 20))

cat("\n=== Script complete ===\n")
