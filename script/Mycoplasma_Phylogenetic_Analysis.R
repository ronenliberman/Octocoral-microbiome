# =============================================================================
# Mycoplasma Phylogenetic Analysis Pipeline
# =============================================================================
# Purpose:  Extract all Mycoplasma ASVs, BLAST against NCBI, build
#           phylogenetic tree with reference sequences, and create 
#           haplotype network across host species
# Author:   Ronen Liberman
# Date:     March 2026
# Manuscript: FEMSEC-25-12-0321 (Revision)
# 
# NOTE: This script follows the same structure as the Endozoicomonas
#       analysis. It assumes you have already loaded and filtered the
#       main dataset (object 'd') from the Endozoicomonas script.
#       If running standalone, uncomment PHASE 0 below.
# =============================================================================

# ---- PACKAGES ----------------------------------------------------------------

library(microeco)
library(file2meco)
library(qiime2R)
library(tidyverse)
library(magrittr)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(ape)
library(pegas)
library(phangorn)
library(msa)
library(httr)
library(xml2)
library(writexl)
library(ggforce)
library(ggrepel)
library(patchwork)
library(vegan)

set.seed(123)

# ---- PATHS -------------------------------------------------------------------

base_dir   <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/16 rRNA analysis/"
rev_dir    <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/"
output_dir <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Output/"
phylo_dir  <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/raw_data/phylo analysis"

myco_output <- file.path(output_dir, "Mycoplasma_phylogeny")
dir.create(myco_output, showWarnings = FALSE, recursive = TRUE)

# ---- COLOR SCHEME (consistent with manuscript PCoA Figure 2) -----------------

species_colors <- c(
  "B. asbestinum"  = "#E41A1C",
  "E. caribaeorum" = "#FF7F00",
  "E. flexuosa"    = "#984EA3",
  "M. muricata"    = "#4DAF4A"
)


# =============================================================================
# PHASE 0: Load data (UNCOMMENT if running standalone, SKIP if 'd' exists)
# =============================================================================

# data <- qiime2meco(
#   file.path(base_dir, "raw data/table-NoEMC.qza"),
#   taxonomy_table = file.path(base_dir, "raw data/taxonomy.qza"),
#   sample_table   = file.path(base_dir, "raw data/metadata_all_edited.csv"),
#   auto_tidy = TRUE
# )
# data$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# d <- clone(data)
# d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE"))
# d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
# d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))
# samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024")
# d$sample_table <- d$sample_table %>% filter(!rownames(.) %in% samples_to_remove)
# species.names <- list("BAST" = "B. asbestinum", "ECAR" = "E. caribaeorum",
#                       "MMUR" = "M. muricata", "EFLE" = "E. flexuosa")
# d$sample_table <- d$sample_table %>% mutate(SPP = recode(SPP, !!!species.names))
# d$tidy_dataset()
# rep_seqs <- readDNAStringSet(
#   file.path(rev_dir, "qiime/rep-seqs-PEE2/dna-sequences-PEE2.fasta")
# )


# =============================================================================
# PHASE 1: Extract all Mycoplasma ASVs
# =============================================================================

cat("\n========== PHASE 1: Extracting Mycoplasma ASVs ==========\n\n")

# --- 1A. Identify Mycoplasma ASVs ---
# Mycoplasma is classified under:
#   Phylum: Firmicutes (formerly Tenericutes in older SILVA versions)
#   Class:  Mollicutes (in SILVA 138.2, may appear as Bacilli)
#   Order:  Mycoplasmatales
#   Family: Mycoplasmataceae
#   Genus:  Mycoplasma

# Search at multiple taxonomic levels to catch all Mycoplasma-related ASVs
myco_genus <- d$tax_table %>%
  filter(grepl("Mycoplasma", Genus, ignore.case = TRUE)) %>% rownames()

myco_family <- d$tax_table %>%
  filter(grepl("Mycoplasmataceae", Family, ignore.case = TRUE)) %>% rownames()

myco_order <- d$tax_table %>%
  filter(grepl("Mycoplasmatales", Order, ignore.case = TRUE)) %>% rownames()

# Also check for Mollicutes class (may capture unclassified Mycoplasma relatives)
mollicutes <- d$tax_table %>%
  filter(grepl("Mollicutes", Class, ignore.case = TRUE)) %>% rownames()

all_myco_asvs <- unique(c(myco_genus, myco_family, myco_order, mollicutes))

cat("Mycoplasma genus-level ASVs:", length(myco_genus), "\n")
cat("Mycoplasmataceae family-level ASVs:", length(myco_family), "\n")
cat("Mycoplasmatales order-level ASVs:", length(myco_order), "\n")
cat("Mollicutes class-level ASVs:", length(mollicutes), "\n")
cat("Combined unique ASVs:", length(all_myco_asvs), "\n\n")

# Check if ASVs exist in OTU table and rep_seqs
cat("In OTU table:", sum(all_myco_asvs %in% rownames(d$otu_table)), "\n")
cat("In rep_seqs:", sum(all_myco_asvs %in% names(rep_seqs)), "\n")

# --- 1B. Build master table ---

myco_otu <- d$otu_table[rownames(d$otu_table) %in% all_myco_asvs, ]

master_list <- lapply(rownames(myco_otu), function(asv_id) {
  counts <- unlist(myco_otu[asv_id, ])
  present_samples <- names(counts[counts > 0])
  if (length(present_samples) == 0) return(NULL)
  
  meta_sub <- d$sample_table[rownames(d$sample_table) %in% present_samples, ]
  host_species <- sort(unique(meta_sub$SPP))
  sites <- sort(unique(meta_sub$Site))
  years <- sort(unique(meta_sub$Year))
  total_reads <- sum(counts[counts > 0])
  n_samples <- length(present_samples)
  
  sample_totals <- colSums(d$otu_table)
  rel_abund <- counts / sample_totals
  mean_rel_abund <- mean(rel_abund[rel_abund > 0]) * 100
  
  tax <- d$tax_table[asv_id, ]
  sequence <- ifelse(asv_id %in% names(rep_seqs),
                     as.character(rep_seqs[[asv_id]]), NA)
  
  # Dominant host
  host_counts <- table(meta_sub$SPP)
  dominant_host <- names(host_counts)[which.max(host_counts)]
  
  data.frame(
    ASV_ID = asv_id,
    Label = paste0("Myco_", which(rownames(myco_otu) == asv_id)),
    Sequence = sequence,
    Seq_Length = nchar(sequence),
    Host_Species = paste(host_species, collapse = "; "),
    N_Host_Species = length(host_species),
    Dominant_Host = dominant_host,
    Sites = paste(sites, collapse = "; "), N_Sites = length(sites),
    Years = paste(years, collapse = "; "), N_Years = length(years),
    Total_Reads = total_reads, N_Samples = n_samples,
    Mean_RelAbund_Pct = round(mean_rel_abund, 4),
    Kingdom = tax$Kingdom, Phylum = tax$Phylum, Class = tax$Class,
    Order = tax$Order, Family = tax$Family, Genus = tax$Genus,
    stringsAsFactors = FALSE
  )
})

myco_master <- bind_rows(master_list) %>% arrange(desc(Total_Reads))
myco_master$Rank <- 1:nrow(myco_master)

# Relabel top ASVs sequentially
myco_master$Label <- paste0("Myco_", myco_master$Rank)

cat("=== MYCOPLASMA MASTER TABLE ===\n")
cat("Total ASVs:", nrow(myco_master), "\n")
cat("ASVs in >1 host:", sum(myco_master$N_Host_Species > 1), "\n")
cat("ASVs in >1 site:", sum(myco_master$N_Sites > 1), "\n\n")

print(myco_master %>%
        select(Label, Host_Species, N_Samples, Mean_RelAbund_Pct, Genus) %>%
        head(15))

write.csv(myco_master, file.path(myco_output, "Mycoplasma_master_table.csv"),
          row.names = FALSE)

# --- 1C. Export FASTA for manual NCBI BLAST ---

myco_seqs <- rep_seqs[names(rep_seqs) %in% myco_master$ASV_ID]
fasta_names <- sapply(names(myco_seqs), function(id) {
  row <- myco_master[myco_master$ASV_ID == id, ]
  paste0(row$Label, "|", id, "|", gsub("; ", "_", row$Host_Species))
})
names(myco_seqs) <- fasta_names
writeXStringSet(myco_seqs, file.path(myco_output, "Mycoplasma_ASVs_for_BLAST.fasta"))

cat("FASTA exported for BLAST:", length(myco_seqs), "sequences\n")
cat("File:", file.path(myco_output, "Mycoplasma_ASVs_for_BLAST.fasta"), "\n")
cat("\n")
cat("=== MANUAL STEP ===\n")
cat("1. Upload this FASTA to NCBI BLAST: https://blast.ncbi.nlm.nih.gov/\n")
cat("2. Database: 16S ribosomal RNA (Bacteria and Archaea)\n")
cat("3. Max target sequences: 10\n")
cat("4. Download results as Hit Table (CSV)\n")
cat("5. Save as 'Myco_blast_results.csv' in:", myco_output, "\n")
cat("6. Then continue with Phase 2 below.\n")

# >>> STOP HERE - go do the manual BLAST, then continue <<<


# =============================================================================
# PHASE 2: Read BLAST Results & Mine Metadata
# =============================================================================

cat("\n========== PHASE 2: Read BLAST Results ==========\n\n")

# --- 2A. Read BLAST CSV ---

blast_web <- read.csv(
  file.path(myco_output, "Myco_blast_results.csv"),
  header = FALSE,
  col.names = c("qseqid", "saccver", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore")
)

blast_local <- blast_web %>%
  mutate(
    Label  = sub("\\|.*", "", qseqid),
    ASV_ID = sub("^[^|]+\\|([^|]+)\\|.*", "\\1", qseqid)
  )

cat("Total BLAST hits:", nrow(blast_local), "\n")
cat("Unique ASVs with hits:", length(unique(blast_local$ASV_ID)), "\n")

# --- 2B. Fetch metadata from NCBI ---

fetch_ncbi_metadata <- function(accession) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                "db=nucleotide&id=", accession, "&rettype=gb&retmode=xml")
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res)) return(list(host = NA, location = NA,
                                source_organism = NA, isolation_source = NA))
  xml_text_content <- content(res, as = "text", encoding = "UTF-8")
  xml_doc <- tryCatch(read_xml(xml_text_content), error = function(e) NULL)
  if (is.null(xml_doc)) return(list(host = NA, location = NA,
                                    source_organism = NA, isolation_source = NA))
  qualifiers <- xml_find_all(xml_doc, ".//GBQualifier")
  host <- NA; location <- NA; source_org <- NA; isolation_source <- NA
  for (q in qualifiers) {
    name  <- xml_text(xml_find_first(q, ".//GBQualifier_name"))
    value <- xml_text(xml_find_first(q, ".//GBQualifier_value"))
    if (!is.na(name)) {
      if (name == "host")             host <- value
      if (name == "country")          location <- value
      if (name == "organism")         source_org <- value
      if (name == "isolation_source") isolation_source <- value
    }
  }
  list(host = host, location = location, source_organism = source_org,
       isolation_source = isolation_source)
}

unique_accessions <- unique(blast_local$saccver)
cat("Fetching metadata for", length(unique_accessions), "accessions...\n")

metadata_list <- list()
for (i in seq_along(unique_accessions)) {
  acc <- unique_accessions[i]
  if (i %% 20 == 0) cat("  Progress:", i, "/", length(unique_accessions), "\n")
  meta <- tryCatch(fetch_ncbi_metadata(acc),
                   error = function(e) list(host = NA, location = NA,
                                            source_organism = NA, isolation_source = NA))
  metadata_list[[acc]] <- data.frame(
    saccver = acc, NCBI_Host = meta$host, NCBI_Location = meta$location,
    NCBI_Organism = meta$source_organism, NCBI_Isolation = meta$isolation_source,
    stringsAsFactors = FALSE
  )
  Sys.sleep(0.4)
}

ncbi_metadata <- bind_rows(metadata_list)
blast_annotated <- blast_local %>% left_join(ncbi_metadata, by = "saccver")

write.csv(blast_annotated, file.path(myco_output, "Mycoplasma_BLAST_annotated.csv"),
          row.names = FALSE)

# --- 2C. Summarize ---

blast_summary <- blast_annotated %>%
  group_by(ASV_ID) %>%
  summarize(
    Best_Hit_Accession = saccver[which.max(bitscore)],
    Best_Hit_Pident    = max(pident),
    Best_Hit_Organism  = NCBI_Organism[which.max(bitscore)],
    Best_Hit_Host      = NCBI_Host[which.max(bitscore)],
    Best_Hit_Location  = NCBI_Location[which.max(bitscore)],
    All_Hosts_Found    = paste(unique(na.omit(NCBI_Host)), collapse = "; "),
    N_Coral_Hits       = sum(grepl(
      "coral|gorgo|octocoral|alcyon|briareum|eunicea|muricea|erythropodium|gorgonia|antillogorgia|plexaura",
      NCBI_Host, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  )

myco_master_annotated <- myco_master %>% left_join(blast_summary, by = "ASV_ID")
write.csv(myco_master_annotated,
          file.path(myco_output, "Mycoplasma_master_BLAST_annotated.csv"),
          row.names = FALSE)

cat("\n=== BLAST SUMMARY ===\n")
cat("Mycoplasma species matched:\n")
blast_annotated %>%
  filter(grepl("Mycoplasma", NCBI_Organism, ignore.case = TRUE)) %>%
  count(NCBI_Organism, sort = TRUE) %>% print()

cat("\nCore ASV best hits:\n")
myco_master_annotated %>%
  head(10) %>%
  select(Label, Best_Hit_Pident, Best_Hit_Organism, Best_Hit_Host) %>%
  print()


# =============================================================================
# PHASE 3: Download Reference Sequences
# =============================================================================

cat("\n========== PHASE 3: Reference Sequences ==========\n\n")

# Mycoplasma type strains relevant to coral/marine systems
# NOTE: Update these accessions based on your BLAST results.
#       The species below are commonly found in coral/octocoral studies.
#       Check your Phase 2 output (NCBI_Organism column) for the specific 
#       species your ASVs match, then find their type strain accessions.

type_strains <- tribble(
  ~Accession, ~Label,                            ~Host,                    ~Location,
  # Common Mycoplasma from coral/marine systems - UPDATE ACCESSIONS
  "AF412966", "M. penetrans HF-2T",              "Human",                  "Reference",
  "NR_025989", "M. iowae 695T",                  "Poultry",                "Reference",
  "NR_113669", "M. todarodis MT1T",              "Squid",                  "Japan",
  "NR_044658", "M. phocicerebrale 1049T",        "Seal",                   "North Sea",
  "NR_113899", "M. insons M1305T",               "Human",                  "Reference"
  # ADD MORE based on your BLAST results:
  # Look at the Best_Hit_Organism column in Mycoplasma_master_BLAST_annotated.csv
  # For each unique Mycoplasma species, find the type strain accession on NCBI
)

# Outgroup: use a non-Mycoplasma Mollicutes or Firmicutes
outgroup <- tribble(
  ~Accession, ~Label,                        ~Host,          ~Location,
  "NR_025949", "Acholeplasma laidlawii T",   "Compost",      "Reference"
)

all_refs <- bind_rows(type_strains, outgroup)
cat("NOTE: Update type_strains with species from your BLAST results!\n")
cat("Reference sequences to download:", nrow(all_refs), "\n\n")

# Download
download_ncbi_fasta <- function(accession) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                "db=nucleotide&id=", accession, "&rettype=fasta&retmode=text")
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  fasta_text <- content(res, as = "text", encoding = "UTF-8")
  lines <- str_split(fasta_text, "\n")[[1]]
  sequence <- paste(lines[-1][lines[-1] != ""], collapse = "")
  if (nchar(sequence) < 100) return(NULL)
  sequence
}

ref_sequences <- list()
for (i in 1:nrow(all_refs)) {
  cat("  Downloading:", all_refs$Accession[i], "-", all_refs$Label[i], "\n")
  seq <- download_ncbi_fasta(all_refs$Accession[i])
  if (!is.null(seq)) {
    ref_sequences[[all_refs$Label[i]]] <- seq
    cat("    ->", nchar(seq), "bp\n")
  }
  Sys.sleep(0.5)
}

cat("Downloaded:", length(ref_sequences), "/", nrow(all_refs), "\n")

# Trim to V4 region
trim_to_v4_fixed <- function(full_seq) {
  f_match <- regexpr("CAGC.GCCGCGGTAA", full_seq, ignore.case = TRUE)
  if (f_match == -1) return(NA)
  r_match <- regexpr("GTAGTCC", full_seq, ignore.case = TRUE)
  if (r_match == -1) return(NA)
  start <- f_match + attr(f_match, "match.length")
  end   <- r_match - 1
  if (end <= start) return(NA)
  if ((end - start) < 200) return(NA)
  if ((end - start) > 300) return(NA)
  substr(full_seq, start, end)
}

ref_v4_list <- list()
for (nm in names(ref_sequences)) {
  v4 <- trim_to_v4_fixed(ref_sequences[[nm]])
  if (!is.na(v4)) {
    ref_v4_list[[nm]] <- v4
    cat("  Trimmed:", nm, "->", nchar(v4), "bp\n")
  } else {
    cat("  WARN: Could not extract V4 from:", nm, "\n")
  }
}
cat("Trimmed successfully:", length(ref_v4_list), "/", length(ref_sequences), "\n")
# =============================================================================
# PHASE 4A: Phylogenetic Tree (NJ + Bootstrap)
# =============================================================================

cat("\n========== PHASE 4A: Phylogenetic Tree ==========\n\n")

# Combine your ASVs + V4-trimmed references
ref_v4_dna <- DNAStringSet(unlist(ref_v4_list))
names(ref_v4_dna) <- paste0("REF|", names(ref_v4_list))

your_seqs <- rep_seqs[names(rep_seqs) %in% myco_master$ASV_ID]
your_labels <- sapply(names(your_seqs), function(id) {
  row <- myco_master[myco_master$ASV_ID == id, ]
  paste0(row$Label, "|", gsub("; ", "+", row$Host_Species))
})
names(your_seqs) <- your_labels

combined_seqs <- c(your_seqs, ref_v4_dna)
writeXStringSet(combined_seqs, file.path(myco_output, "Mycoplasma_combined_for_tree.fasta"))
cat("Combined sequences:", length(combined_seqs), "\n")

# Align
cat("Running alignment...\n")
aligned <- msa(combined_seqs, method = "Muscle")
aligned_dnabin <- as.DNAbin(aligned)
cat("Alignment complete.\n")

# Save alignment (for Galaxy/RAxML if desired)
write.FASTA(aligned_dnabin, file.path(myco_output, "Mycoplasma_alignment.fasta"))

# Also save clean version for RAxML (no special characters in names)
clean_names <- labels(aligned_dnabin)
clean_names <- gsub("\\|", "_", clean_names)
clean_names <- gsub("\\+", "_", clean_names)
clean_names <- gsub(" ", "_", clean_names)
clean_names <- gsub("\\.", "", clean_names)
clean_names <- gsub(",", "", clean_names)
clean_names <- make.unique(clean_names, sep = "_")
aligned_clean <- aligned_dnabin
rownames(aligned_clean) <- clean_names
write.FASTA(aligned_clean, file.path(myco_output, "Mycoplasma_alignment_CLEAN.fasta"))

# Name mapping for converting back
name_map <- data.frame(Original = labels(aligned_dnabin), Clean = clean_names,
                       stringsAsFactors = FALSE)
write.csv(name_map, file.path(myco_output, "Mycoplasma_name_mapping.csv"), row.names = FALSE)

# Build NJ tree
phy_dat <- as.phyDat(aligned_dnabin)
dm <- dist.ml(phy_dat, model = "JC69")
nj_tree <- NJ(dm)

# Root with outgroup (adjust label if needed)
outgroup_label <- grep("Acholeplasma", nj_tree$tip.label, value = TRUE)
if (length(outgroup_label) > 0) {
  nj_tree <- root(nj_tree, outgroup = outgroup_label[1], resolve.root = TRUE)
}

# Bootstrap
cat("Running 1000 NJ bootstrap replicates...\n")
set.seed(123)
boot_trees <- boot.phylo(nj_tree, aligned_dnabin,
                         function(x) NJ(dist.dna(x, model = "JC69")),
                         B = 1000, quiet = FALSE)

tree_bs <- nj_tree
tree_bs$node.label <- boot_trees / 10  # Convert to percentages

write.tree(tree_bs, file.path(myco_output, "Mycoplasma_NJ_bootstrap.nwk"))
cat("Tree saved.\n")


# =============================================================================
# # PHASE 4A: MYCOPLASMA IQ-TREE Visualization
# =============================================================================
options(ignore.negative.edge = TRUE)
# Read tree
myco_iqtree <- read.tree(file.path(phylo_dir, "Myco_iqtree.treefile"))

# Root with outgroup
outgroup_tip <- grep("Acholeplasma", myco_iqtree$tip.label, value = TRUE)
cat("Outgroup:", outgroup_tip, "\n")
myco_iqtree <- root(myco_iqtree, outgroup = outgroup_tip, resolve.root = TRUE)

# Parse UFboot values
myco_iqtree$node.label <- sapply(myco_iqtree$node.label, function(x) {
  if (is.na(x) || x == "") return(NA)
  parts <- strsplit(as.character(x), "/")[[1]]
  if (length(parts) == 2) return(as.numeric(parts[2]))
  return(NA)
})

tree_bs_myco <- myco_iqtree

# Recreate Dominant_Host_Filtered if not present
if (!"Dominant_Host_Filtered" %in% colnames(myco_master)) {
  sample_totals <- colSums(d$otu_table)
  myco_master$Dominant_Host_Filtered <- sapply(myco_master$ASV_ID, function(asv_id) {
    counts <- unlist(d$otu_table[asv_id, ])
    present_samples <- names(counts[counts > 0])
    host_counts_filtered <- c()
    for (sp in c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")) {
      sp_samples <- rownames(d$sample_table[d$sample_table$SPP == sp, ])
      sp_present <- intersect(present_samples, sp_samples)
      count <- 0
      for (samp in sp_present) {
        total_myco <- sum(unlist(myco_otu[, samp]))
        total_reads <- sample_totals[samp]
        if ((total_myco / total_reads) > 0.01) count <- count + 1
      }
      host_counts_filtered[sp] <- count
    }
    if (max(host_counts_filtered) == 0) return("Low abundance")
    names(host_counts_filtered)[which.max(host_counts_filtered)]
  })
}

# Build tip_data
tip_data_myco <- data.frame(
  label = tree_bs_myco$tip.label,
  stringsAsFactors = FALSE
)
tip_data_myco$Myco_Label <- sub("\\|.*", "", tip_data_myco$label)
tip_data_myco$Source <- ifelse(grepl("^REF", tip_data_myco$label), "Type strain", "This study")

top_n_labels <- myco_master$Label[1:min(10, nrow(myco_master))]
tip_data_myco$Short_Label <- ifelse(
  tip_data_myco$Myco_Label %in% top_n_labels,
  tip_data_myco$Myco_Label,
  ifelse(grepl("^REF", tip_data_myco$label),
         gsub("_", " ", sub("^REF\\|", "", tip_data_myco$label)),
         "")
)

# Host_Group using 1% filtered dominant host
tip_data_myco$Host_Group <- NA_character_
for (i in 1:nrow(tip_data_myco)) {
  lbl <- tip_data_myco$Myco_Label[i]
  if (grepl("^REF", lbl)) {
    tip_data_myco$Host_Group[i] <- "Reference"
  } else {
    idx <- which(myco_master$Label == lbl)
    if (length(idx) > 0) {
      tip_data_myco$Host_Group[i] <- myco_master$Dominant_Host_Filtered[idx[1]]
    } else {
      tip_data_myco$Host_Group[i] <- "Unknown"
    }
  }
}

cat("Myco Host_Group values:\n")
print(table(tip_data_myco$Host_Group))

# Colors
tree_colors_myco <- c(
  "B. asbestinum" = "#E41A1C", "E. caribaeorum" = "#FF7F00",
  "E. flexuosa" = "#984EA3", "M. muricata" = "#4DAF4A",
  "Reference" = "black", "Unknown" = "grey50",
  "Low abundance" = "grey80"
)

# Plot
p_tree_myco <- ggtree(tree_bs_myco, layout = "rectangular", size = 0.5, color = "grey30") %<+%
  tip_data_myco +
  geom_tippoint(aes(color = Host_Group, shape = Source), size = 3, stroke = 0.5) +
  geom_tiplab(aes(label = Short_Label), size = 3.5, offset = 0.002, fontface = "bold") +
  geom_nodelab(aes(label = ifelse(!is.na(as.numeric(label)) & as.numeric(label) >= 70,
                                  round(as.numeric(label)), "")),
               size = 2.5, hjust = -0.2, vjust = -0.5, color = "black") +
  scale_color_manual(
    values = tree_colors_myco, name = "Dominant Host",
    breaks = c("B. asbestinum", "E. flexuosa", "M. muricata", "Low abundance", "Reference"),
    labels = c(expression(italic("B. asbestinum")),
               expression(italic("E. flexuosa")),
               expression(italic("M. muricata")),
               "Low abundance (<1%)", "Reference type strain")
  ) +
  scale_shape_manual(values = c("This study" = 16, "Type strain" = 17), name = "Source") +
  geom_treescale(x = 0, y = -2, fontsize = 3.5, linesize = 0.8, offset = 1) +
  theme_tree() +
  theme(
    legend.position = c(0.85, 0.35),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    plot.margin = margin(10, 30, 10, 10)
  ) +
  annotate("text", x = -Inf, y = Inf, label = "A",
           size = 8, fontface = "bold", hjust = -0.5, vjust = 1.5)

ggsave(file.path(myco_output, "Mycoplasma_IQtree.png"),
       p_tree_myco, width = 12, height = 14, dpi = 300, bg = "white")

cat("Mycoplasma IQ-TREE figure saved.\n")
print(p_tree_myco)

# =============================================================================
# PHASE 4B: Haplotype Network
# =============================================================================

# Mycoplasma haplotype network - run all at once
keep_asvs <- myco_master %>% filter(N_Samples >= 3) %>% pull(Label)
filtered_seqs <- myco_seqs_for_haplo[names(myco_seqs_for_haplo) %in% keep_asvs]
cat("Filtered ASVs for network:", length(filtered_seqs), "\n")

filtered_aligned <- msa(filtered_seqs, method = "Muscle")
filtered_bin <- as.DNAbin(filtered_aligned)

hap_myco <- haplotype(filtered_bin)
net_myco <- haploNet(hap_myco)
hap_index_myco <- attr(hap_myco, "index")

host_species_levels <- c("B. asbestinum", "E. caribaeorum",
                         "E. flexuosa", "M. muricata")

# 1% filtered frequency matrix
sample_totals <- colSums(d$otu_table)
freq_matrix_myco_filtered <- matrix(0, nrow = length(hap_index_myco), ncol = 4,
                                    dimnames = list(NULL, host_species_levels))

for (h in seq_along(hap_index_myco)) {
  for (idx in hap_index_myco[[h]]) {
    asv_label <- labels(filtered_bin)[idx]
    row <- myco_master[myco_master$Label == asv_label, ]
    if (nrow(row) == 0) next
    asv_id <- row$ASV_ID
    counts <- unlist(d$otu_table[asv_id, ])
    present_samples <- names(counts[counts > 0])
    for (sp in host_species_levels) {
      sp_samples <- rownames(d$sample_table[d$sample_table$SPP == sp, ])
      sp_present <- intersect(present_samples, sp_samples)
      for (samp in sp_present) {
        total_myco <- sum(unlist(myco_otu[, samp]))
        total_reads <- sample_totals[samp]
        if ((total_myco / total_reads) > 0.01) {
          freq_matrix_myco_filtered[h, sp] <- freq_matrix_myco_filtered[h, sp] + 1
        }
      }
    }
  }
}

hap_sizes_myco_filtered <- rowSums(freq_matrix_myco_filtered)

# Labels for top ASVs
top_n_labels <- myco_master$Label[1:min(10, nrow(myco_master))]
hap_label_vec_myco <- rep("", length(hap_index_myco))
for (h in seq_along(hap_index_myco)) {
  for (idx in hap_index_myco[[h]]) {
    lbl <- labels(filtered_bin)[idx]
    if (lbl %in% top_n_labels) hap_label_vec_myco[h] <- lbl
  }
}

cat("Haplotypes:", length(hap_sizes_myco_filtered), "\n")
cat("Per species:\n")
print(colSums(freq_matrix_myco_filtered))

# Extract coordinates
png(tempfile())
coords <- plot(net_myco, size = sqrt(hap_sizes_myco_filtered + 1) * 1.5,
               pie = freq_matrix_myco_filtered,
               bg = species_colors, labels = FALSE, scale.ratio = 0.3,
               show.mutation = 1, threshold = 0)
dev.off()

node_x <- coords$xx
node_y <- coords$yy
node_size <- coords$size
pie_data <- coords$pie
links <- coords$link
alter_links <- coords$alter.links

edge_df <- data.frame(
  x = node_x[links[, 1]], y = node_y[links[, 1]],
  xend = node_x[links[, 2]], yend = node_y[links[, 2]]
)

mut_df <- data.frame(
  x = alter_links[, 1], y = alter_links[, 2],
  xend = alter_links[, 3], yend = alter_links[, 4]
)

pie_list <- list()
for (i in 1:nrow(pie_data)) {
  props <- pie_data[i, ]
  total <- sum(props)
  if (total == 0) next
  fracs <- props / total
  end_angles <- cumsum(fracs) * 2 * pi
  start_angles <- c(0, end_angles[-length(end_angles)])
  for (j in 1:length(fracs)) {
    if (fracs[j] == 0) next
    pie_list[[length(pie_list) + 1]] <- data.frame(
      x0 = node_x[i], y0 = node_y[i], r = node_size[i] * 0.55,
      start = start_angles[j], end = end_angles[j],
      species = colnames(pie_data)[j], stringsAsFactors = FALSE
    )
  }
}
pie_df <- bind_rows(pie_list)

label_df <- data.frame(x = node_x, y = node_y, label = hap_label_vec_myco,
                       stringsAsFactors = FALSE) %>% filter(label != "")
  
cat("Haplotype network data ready. Now run Figure 6 code.\n")
  
  p_haplo_filtered <- ggplot() +
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey40", linewidth = 0.5) +
    geom_segment(data = mut_df, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 0.6) +
    geom_arc_bar(data = pie_df,
                 aes(x0 = x0, y0 = y0, r0 = 0, r = r,
                     start = start, end = end, fill = species),
                 color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = species_colors, name = "Host Species",
      labels = c(expression(italic("B. asbestinum")),
                 expression(italic("E. caribaeorum")),
                 expression(italic("E. flexuosa")),
                 expression(italic("M. muricata")))
    ) +
    geom_label_repel(data = label_df, aes(x = x, y = y, label = label),
                     size = 4, fontface = "bold", fill = "white", alpha = 0.85,
                     box.padding = 0.8, segment.color = "grey30",
                     min.segment.length = 0, max.overlaps = 20) +
    annotate("text", x = -Inf, y = Inf, label = "B",
             size = 8, fontface = "bold", hjust = -0.5, vjust = 1.5) +
    coord_equal() + theme_void() +
    theme(
      legend.position = c(0.88, 0.85),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.6, "cm"),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  ggsave(file.path(myco_output, "Mycoplasma_haplotype_network_filtered.png"),
         p_haplo_filtered, width = 12, height = 10, bg= "white", dpi = 300)
  
  cat("Filtered haplotype network saved.\n")
  print(p_haplo_filtered)
  
  
  # =============================================================================
  # FIGURE 6: Combined Mycoplasma Tree + Haplotype Network
  # =============================================================================
  
  # Tree panel - no individual tag, legend at bottom
  p_tree_fig6 <- ggtree(tree_bs_myco, layout = "rectangular", size = 0.5, color = "grey30") %<+%
    tip_data_myco +
    geom_tippoint(aes(color = Host_Group, shape = Source), size = 3, stroke = 0.5) +
    geom_tiplab(aes(label = Short_Label), size = 4.5, offset = 0.005, fontface = "bold") +
    geom_nodelab(aes(label = ifelse(!is.na(as.numeric(label)) & as.numeric(label) >= 70,
                                    round(as.numeric(label)), "")),
                 size = 2.5, hjust = -0.2, vjust = -0.5, color = "black") +
    scale_color_manual(
      values = tree_colors_myco, name = "Dominant Host",
      breaks = c("B. asbestinum", "E. flexuosa", "M. muricata", "Low abundance", "Reference"),
      labels = c(expression(italic("B. asbestinum")),
                 expression(italic("E. flexuosa")),
                 expression(italic("M. muricata")),
                 "Low abundance (<1%)", "Reference type strain")
    ) +
    scale_shape_manual(values = c("This study" = 16, "Type strain" = 17), name = "Source") +
    geom_treescale(x = max(fortify(tree_bs_myco)$x) * 0.8, y = -2,
                   fontsize = 5, linesize = 1.5, offset = 1) +
    theme_tree() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.6, "cm"),
      plot.margin = margin(5, 0, 5, 5)
    ) +
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 1))
  
  # Haplotype network - no legend (tree legend covers both)
  p_haplo_fig6 <- ggplot() +
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey40", linewidth = 0.5) +
    geom_segment(data = mut_df, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 0.6) +
    geom_arc_bar(data = pie_df,
                 aes(x0 = x0, y0 = y0, r0 = 0, r = r,
                     start = start, end = end, fill = species),
                 color = "black", linewidth = 0.3) +
    scale_fill_manual(values = species_colors, guide = "none") +
    geom_label_repel(data = label_df, aes(x = x, y = y, label = label),
                     size = 4.5, fontface = "bold", fill = "white", alpha = 0.85,
                     box.padding = 0.8, segment.color = "grey30",
                     min.segment.length = 0, max.overlaps = 20) +
    coord_equal() + theme_void() +
    theme(plot.margin = margin(5, 5, 5, 0))
  
  # Combine
  fig6 <- p_tree_fig6 + p_haplo_fig6 +
    plot_layout(widths = c(1, 1.2)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  ggsave(file.path(myco_output, "Figure6_Mycoplasma_combined.png"),
         fig6, width = 14, height = 14, dpi = 300, bg = "white")
  ggsave(file.path(myco_output, "Figure6_Mycoplasma_combined.pdf"),
         fig6, width = 20, height = 14, bg = "white")
  
  cat("Figure 6 saved.\n")
  print(fig6)