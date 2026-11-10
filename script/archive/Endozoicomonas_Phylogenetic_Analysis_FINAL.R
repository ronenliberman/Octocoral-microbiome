# =============================================================================
# Endozoicomonas Phylogenetic Analysis Pipeline - FINAL CLEAN VERSION
# =============================================================================
# Purpose:  Extract all Endozoicomonas ASVs, BLAST against NCBI, build
#           ML phylogenetic tree (IQ-TREE via terminal), and create 
#           haplotype network across host species
# Author:   Ronen Liberman
# Date:     March 2026
# Manuscript: FEMSEC-25-12-0321 (Revision)
#
# WORKFLOW OVERVIEW:
#   Phase 1: Extract Endozoicomonas ASVs → master table + FASTA for BLAST
#   Phase 2: Read NCBI BLAST results → annotate with host/location metadata
#   Phase 3: Download reference type strain sequences → V4 trim → combine
#   --- TERMINAL STEP: Run MAFFT + IQ-TREE on Mac ---
#   Phase 4A: Visualize IQ-TREE ML tree with ggtree
#   Phase 4B: Haplotype network with ggplot2 + ggforce
#   Phase 5: Site-specific PERMANOVA tests
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

endo_output <- file.path(output_dir, "Endozoicomonas_phylogeny")
dir.create(endo_output, showWarnings = FALSE, recursive = TRUE)

# ---- COLOR SCHEME (consistent with manuscript PCoA Figure 2) -----------------

species_colors <- c(
  "B. asbestinum"  = "#E41A1C",
  "E. caribaeorum" = "#FF7F00",
  "E. flexuosa"    = "#984EA3",
  "M. muricata"    = "#4DAF4A"
)


# =============================================================================
# PHASE 1: Extract all Endozoicomonas ASVs
# =============================================================================

cat("\n========== PHASE 1: Extracting Endozoicomonas ASVs ==========\n\n")

# --- 1A. Load and filter data ---

data <- qiime2meco(
  file.path(base_dir, "raw data/table-NoEMC.qza"),
  taxonomy_table = file.path(base_dir, "raw data/taxonomy.qza"),
  sample_table   = file.path(base_dir, "raw data/metadata_all_edited.csv"),
  auto_tidy = TRUE
)

data$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

d <- clone(data)
d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE"))
d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

# Remove unreliable samples
samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024", "EFLE_DC3_3_2024")
d$sample_table <- d$sample_table %>%
  filter(!rownames(.) %in% samples_to_remove)

species.names <- list(
  "BAST" = "B. asbestinum", "ECAR" = "E. caribaeorum",
  "MMUR" = "M. muricata",  "EFLE" = "E. flexuosa"
)
d$sample_table <- d$sample_table %>% mutate(SPP = recode(SPP, !!!species.names))
d$tidy_dataset()

cat("Filtered dataset:", nrow(d$sample_table), "samples,", 
    nrow(d$otu_table), "ASVs\n")

# --- 1B. Read representative sequences ---

rep_seqs <- readDNAStringSet(
  file.path(rev_dir, "qiime/rep-seqs-PEE2/dna-sequences-PEE2.fasta")
)
cat("Representative sequences loaded:", length(rep_seqs), "\n")

# --- 1C. Identify all Endozoicomonas/Endozoicomonadaceae ASVs ---

endo_asvs <- d$tax_table %>%
  filter(grepl("Endozoicomonas", Genus, ignore.case = TRUE)) %>% rownames()
endoFamily_asvs <- d$tax_table %>%
  filter(grepl("Endozoicomonadaceae", Family, ignore.case = TRUE)) %>% rownames()
all_endo_asvs <- unique(c(endo_asvs, endoFamily_asvs))

cat("Total Endozoicomonas/Endozoicomonadaceae ASVs:", length(all_endo_asvs), "\n")

# --- 1D. Build master table ---

endo_otu <- d$otu_table[rownames(d$otu_table) %in% all_endo_asvs, ]

core_endo_ids <- c(
  "5656aae49e1f57ea6217094b1969e8aa", "238ab607aedd1bd3b1b75b03f61331c3",
  "47ce1658b2ee5601bd633f5f2d014adf", "2f40b4eaa74d3e96797cbb513a5f43bd",
  "a80ddcda6c71332a7a6269bdce4dee87", "69d03d98367f96a5ec36f2c90e977fad"
)
endo_labels_map <- c(
  "5656aae49e1f57ea6217094b1969e8aa" = "Endo_4",
  "238ab607aedd1bd3b1b75b03f61331c3" = "Endo_1",
  "47ce1658b2ee5601bd633f5f2d014adf" = "Endo_3",
  "2f40b4eaa74d3e96797cbb513a5f43bd" = "Endo_2",
  "a80ddcda6c71332a7a6269bdce4dee87" = "Endo_6",
  "69d03d98367f96a5ec36f2c90e977fad" = "Endo_5"
)

master_list <- lapply(rownames(endo_otu), function(asv_id) {
  counts <- unlist(endo_otu[asv_id, ])
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
  
  core_status <- ifelse(asv_id %in% core_endo_ids, "Yes", "No")
  label <- ifelse(asv_id %in% names(endo_labels_map), endo_labels_map[asv_id],
                  paste0("Endo_rare_", which(rownames(endo_otu) == asv_id)))
  
  # Dominant host (by sample count)
  host_counts <- table(meta_sub$SPP)
  dominant_host <- names(host_counts)[which.max(host_counts)]
  
  data.frame(
    ASV_ID = asv_id, Label = label, Sequence = sequence,
    Seq_Length = nchar(sequence),
    Host_Species = paste(host_species, collapse = "; "),
    N_Host_Species = length(host_species),
    Dominant_Host = dominant_host,
    Sites = paste(sites, collapse = "; "), N_Sites = length(sites),
    Years = paste(years, collapse = "; "), N_Years = length(years),
    Total_Reads = total_reads, N_Samples = n_samples,
    Mean_RelAbund_Pct = round(mean_rel_abund, 4),
    Core_90pct = core_status,
    Kingdom = tax$Kingdom, Phylum = tax$Phylum, Class = tax$Class,
    Order = tax$Order, Family = tax$Family, Genus = tax$Genus,
    stringsAsFactors = FALSE
  )
})

endo_master <- bind_rows(master_list) %>% arrange(desc(Total_Reads))
endo_master$Rank <- 1:nrow(endo_master)

cat("=== ENDOZOICOMONAS MASTER TABLE ===\n")
cat("Total ASVs:", nrow(endo_master), "\n")
cat("Core ASVs:", sum(endo_master$Core_90pct == "Yes"), "\n")
cat("ASVs in >1 host:", sum(endo_master$N_Host_Species > 1), "\n")
cat("ASVs in >1 site:", sum(endo_master$N_Sites > 1), "\n\n")

write.csv(endo_master, file.path(endo_output, "Endozoicomonas_master_table.csv"),
          row.names = FALSE)

# --- 1E. Export FASTA for BLAST ---

endo_seqs <- rep_seqs[names(rep_seqs) %in% endo_master$ASV_ID]
new_names <- sapply(names(endo_seqs), function(id) {
  row <- endo_master[endo_master$ASV_ID == id, ]
  paste0(row$Label, "|", id, "|", gsub("; ", "_", row$Host_Species))
})
names(endo_seqs) <- new_names
writeXStringSet(endo_seqs, file.path(endo_output, "Endozoicomonas_ASVs_for_BLAST.fasta"))
cat("FASTA exported for BLAST:", length(endo_seqs), "sequences\n")


# =============================================================================
# PHASE 2: Read NCBI BLAST Results & Mine Metadata
# =============================================================================
# MANUAL STEP: Upload Endozoicomonas_ASVs_for_BLAST.fasta to NCBI BLAST
#   https://blast.ncbi.nlm.nih.gov/ → Nucleotide BLAST
#   Database: 16S ribosomal RNA → Max target sequences: 10
#   Download Hit Table (CSV) → save as "Endo_blast_results.csv" in endo_output

cat("\n========== PHASE 2: Read BLAST Results ==========\n\n")

blast_web <- read.csv(
  file.path(endo_output, "Endo_blast_results.csv"),
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

# --- 2B. Fetch host/location metadata from NCBI ---

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
write.csv(blast_annotated, file.path(endo_output, "Endozoicomonas_BLAST_annotated.csv"),
          row.names = FALSE)

# --- 2C. Summarize per ASV ---

blast_summary <- blast_annotated %>%
  group_by(ASV_ID) %>%
  summarize(
    Best_Hit_Accession = saccver[which.max(bitscore)],
    Best_Hit_Pident    = max(pident),
    Best_Hit_Organism  = NCBI_Organism[which.max(bitscore)],
    Best_Hit_Host      = NCBI_Host[which.max(bitscore)],
    Best_Hit_Location  = NCBI_Location[which.max(bitscore)],
    All_Hosts_Found    = paste(unique(na.omit(NCBI_Host)), collapse = "; "),
    N_Octocoral_Hits   = sum(grepl(
      "coral|gorgo|octocoral|alcyon|briareum|eunicea|muricea|erythropodium|gorgonia|antillogorgia|plexaura|lobophytum|sinularia|eunicella|paramuricea|sarcophyton",
      NCBI_Host, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  )

endo_master_annotated <- endo_master %>% left_join(blast_summary, by = "ASV_ID")
write.csv(endo_master_annotated,
          file.path(endo_output, "Endozoicomonas_master_BLAST_annotated.csv"),
          row.names = FALSE)

cat("=== BLAST SUMMARY ===\n")
cat("All Endozoicomonas species matched:\n")
blast_annotated %>%
  filter(grepl("Endozoicomonas", NCBI_Organism)) %>%
  count(NCBI_Organism, sort = TRUE) %>% print()

cat("\nCore ASV best BLAST hits:\n")
endo_master_annotated %>%
  filter(Core_90pct == "Yes") %>%
  dplyr::select(Label, Host_Species, Best_Hit_Pident, Best_Hit_Organism,
                Best_Hit_Host, N_Octocoral_Hits) %>%
  print(width = 120)


# =============================================================================
# PHASE 3: Download Reference Type Strain Sequences
# =============================================================================

cat("\n========== PHASE 3: Reference Sequences ==========\n\n")

type_strains <- tribble(
  ~Accession, ~Label,                          ~Host,                        ~Location,
  "JX488684", "E. euniceicola EF212T",         "Eunicea fusca",              "Florida, USA",
  "JX488685", "E. gorgoniicola PS125T",        "Plexaura sp.",               "Florida, USA",
  "AB196667", "E. elysicola MKT110T",          "Sea slug Elysia ornata",     "Japan",
  "AB695088", "E. numazuensis HC50T",          "Sponge",                     "Japan",
  "FJ347758", "E. montiporae CL-33T",         "Montipora aequituberculata", "Taiwan",
  "KC878324", "E. atrinae WP70T",             "Atrina pectinata",           "South Korea",
  "KY770540", "E. acroporae Acr-14T",         "Acropora sp.",               "Taiwan",
  "KY770538", "E. ascidiicola AVMART05T",     "Ascidiella sp.",             "Spain",
  "MG456764", "E. lisbonensis 4ANDT",         "Ascidiella sp.",             "Portugal"
)

outgroup <- tribble(
  ~Accession, ~Label,                  ~Host,              ~Location,
  "EU599216", "Kistimonas asteriae T", "Asterias starfish", "South Korea"
)

all_refs <- bind_rows(type_strains, outgroup)

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

cat("Downloading", nrow(all_refs), "reference sequences...\n")
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

# Trim to V4 region (515F-806R)
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

# --- Combine your ASVs + references → export for MAFFT/IQ-TREE ---

ref_v4_dna <- DNAStringSet(unlist(ref_v4_list))
names(ref_v4_dna) <- paste0("REF|", names(ref_v4_list))

your_seqs <- rep_seqs[names(rep_seqs) %in% endo_master$ASV_ID]
your_labels <- sapply(names(your_seqs), function(id) {
  row <- endo_master[endo_master$ASV_ID == id, ]
  paste0(row$Label, "|", gsub("; ", "+", row$Host_Species))
})
names(your_seqs) <- your_labels

combined_seqs <- c(your_seqs, ref_v4_dna)
writeXStringSet(combined_seqs, file.path(endo_output, "Endozoicomonas_combined_for_tree.fasta"))
cat("Combined sequences exported:", length(combined_seqs), "\n")
cat("\n")
cat("=== TERMINAL STEP ===\n")
cat("Transfer Endozoicomonas_combined_for_tree.fasta to Mac, then run:\n")
cat('  mafft --auto "Endozoicomonas_combined_for_tree.fasta" > "Endo_mafft_aligned.fasta"\n')
cat('  iqtree -s "Endo_mafft_aligned.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre "Endo_iqtree"\n')
cat("Transfer Endo_iqtree.treefile back to PC.\n")
cat("Best-fit model: TPM3+I+R3 (check Endo_iqtree.iqtree for confirmation)\n")


# =============================================================================
# PHASE 4A: IQ-TREE ML Phylogenetic Tree Visualization
# =============================================================================

cat("\n========== PHASE 4A: IQ-TREE Visualization ==========\n\n")

options(ignore.negative.edge = TRUE)

# Read IQ-TREE result
tree_bs <- read.tree(file.path(phylo_dir, "Endo_iqtree.treefile"))

# Root with outgroup
outgroup_tip <- grep("Kistimonas", tree_bs$tip.label, value = TRUE)
tree_bs <- root(tree_bs, outgroup = outgroup_tip, resolve.root = TRUE)

# Parse UFboot values from node labels (format: "SH-aLRT/UFboot")
tree_bs$node.label <- sapply(tree_bs$node.label, function(x) {
  if (is.na(x) || x == "") return(NA)
  parts <- strsplit(as.character(x), "/")[[1]]
  if (length(parts) == 2) return(as.numeric(parts[2]))
  return(NA)
})

# Build tip_data
tip_data <- data.frame(label = tree_bs$tip.label, stringsAsFactors = FALSE)
tip_data$Endo_Label <- sub("\\|.*", "", tip_data$label)
tip_data$Source <- ifelse(grepl("^REF", tip_data$label), "Type strain", "This study")
tip_data$Short_Label <- ifelse(
  grepl("^Endo_[1-6]\\|", tip_data$label),
  sub("\\|.*", "", tip_data$label),
  ifelse(grepl("^REF", tip_data$label),
         gsub("_", " ", sub("^REF\\|", "", tip_data$label)),
         "")
)

# Assign Host_Group by matching to endo_master
tip_data$Host_Group <- NA_character_
for (i in 1:nrow(tip_data)) {
  lbl <- tip_data$Endo_Label[i]
  if (grepl("^REF", lbl)) {
    tip_data$Host_Group[i] <- "Reference"
  } else {
    idx <- which(endo_master$Label == lbl)
    if (length(idx) > 0) {
      tip_data$Host_Group[i] <- endo_master$Dominant_Host[idx[1]]
    } else {
      tip_data$Host_Group[i] <- "Unknown"
    }
  }
}

cat("Host_Group distribution:\n")
print(table(tip_data$Host_Group))

# Colors
tree_colors <- c(
  "B. asbestinum" = "#E41A1C", "E. caribaeorum" = "#FF7F00",
  "E. flexuosa" = "#984EA3", "M. muricata" = "#4DAF4A",
  "Reference" = "black", "Unknown" = "grey50"
)

# Plot
p_tree <- ggtree(tree_bs, layout = "rectangular", size = 0.5, color = "grey30") %<+%
  tip_data +
  geom_tippoint(aes(color = Host_Group, shape = Source), size = 3, stroke = 0.5) +
  geom_tiplab(aes(label = Short_Label), size = 3.5, offset = 0.002, fontface = "bold") +
  geom_nodelab(aes(label = ifelse(!is.na(as.numeric(label)) & as.numeric(label) >= 70,
                                  round(as.numeric(label)), "")),
               size = 2.5, hjust = -0.2, vjust = -0.5, color = "black") +
  scale_color_manual(
    values = tree_colors, name = "Dominant Host",
    breaks = c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata", "Reference"),
    labels = c(expression(italic("B. asbestinum")), expression(italic("E. caribaeorum")),
               expression(italic("E. flexuosa")), expression(italic("M. muricata")),
               "Reference type strain")
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

ggsave(file.path(endo_output, "Figure5A_Endo_IQtree.png"),
       p_tree, width = 12, height = 14, dpi = 300, bg = "white")
cat("Figure 5A saved.\n")


# =============================================================================
# PHASE 4B: Haplotype Network (ggplot2 + ggforce)
# =============================================================================

cat("\n========== PHASE 4B: Haplotype Network ==========\n\n")

# Prepare sequences with labels
endo_seqs_for_haplo <- rep_seqs[names(rep_seqs) %in% endo_master$ASV_ID]
haplo_names <- sapply(names(endo_seqs_for_haplo), function(id) {
  endo_master$Label[endo_master$ASV_ID == id]
})
names(endo_seqs_for_haplo) <- haplo_names

# Filter to ASVs with >= 5 samples
keep_asvs <- endo_master %>% filter(N_Samples >= 5) %>% pull(Label)
filtered_seqs <- endo_seqs_for_haplo[names(endo_seqs_for_haplo) %in% keep_asvs]
cat("Filtered ASVs for network:", length(filtered_seqs), "\n")

# Align
filtered_aligned <- msa(filtered_seqs, method = "Muscle")
filtered_bin <- as.DNAbin(filtered_aligned)

# Build haplotype network
hap2 <- haplotype(filtered_bin)
net2 <- haploNet(hap2)
hap_index2 <- attr(hap2, "index")

# Build frequency matrix by host species
host_species_levels <- c("B. asbestinum", "E. caribaeorum",
                         "E. flexuosa", "M. muricata")
freq_matrix2 <- matrix(0, nrow = length(hap_index2), ncol = 4,
                       dimnames = list(NULL, host_species_levels))

for (h in seq_along(hap_index2)) {
  for (idx in hap_index2[[h]]) {
    asv_label <- labels(filtered_bin)[idx]
    row <- endo_master[endo_master$Label == asv_label, ]
    if (nrow(row) == 0) next
    asv_id <- row$ASV_ID
    counts <- unlist(d$otu_table[asv_id, ])
    present_samples <- names(counts[counts > 0])
    meta_sub <- d$sample_table[rownames(d$sample_table) %in% present_samples, ]
    for (sp in host_species_levels) {
      freq_matrix2[h, sp] <- freq_matrix2[h, sp] + sum(meta_sub$SPP == sp)
    }
  }
}

hap_sizes2 <- rowSums(freq_matrix2)

# Label core ASVs
hap_label_vec2 <- rep("", length(hap_index2))
for (h in seq_along(hap_index2)) {
  for (idx in hap_index2[[h]]) {
    lbl <- labels(filtered_bin)[idx]
    if (grepl("^Endo_[1-6]$", lbl)) hap_label_vec2[h] <- lbl
  }
}

cat("Haplotypes:", length(hap_sizes2), "\n")
cat("Core labels assigned:", sum(hap_label_vec2 != ""), "\n")

# --- Extract coordinates and build ggplot2 haplotype network ---

png(tempfile())
coords <- plot(net2, size = sqrt(hap_sizes2) * 1.5, pie = freq_matrix2,
               bg = species_colors, labels = FALSE, scale.ratio = 0.3,
               show.mutation = 1, threshold = 0)
dev.off()

node_x <- coords$xx
node_y <- coords$yy
node_size <- coords$size
pie_data <- coords$pie
links <- coords$link
alter_links <- coords$alter.links

# Edge dataframe
edge_df <- data.frame(
  x = node_x[links[, 1]], y = node_y[links[, 1]],
  xend = node_x[links[, 2]], yend = node_y[links[, 2]]
)

# Mutation tick marks
mut_df <- data.frame(
  x = alter_links[, 1], y = alter_links[, 2],
  xend = alter_links[, 3], yend = alter_links[, 4]
)

# Pie chart wedges
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

# Label dataframe
label_df <- data.frame(
  x = node_x, y = node_y, label = hap_label_vec2,
  stringsAsFactors = FALSE
) %>% filter(label != "")

# Build ggplot
p_haplo <- ggplot() +
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
    labels = c(expression(italic("B. asbestinum")), expression(italic("E. caribaeorum")),
               expression(italic("E. flexuosa")), expression(italic("M. muricata")))
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

ggsave(file.path(endo_output, "Figure5B_Endo_haplotype_network.png"),
       p_haplo, width = 12, height = 10, dpi = 300, bg = "white")
cat("Figure 5B saved.\n")


# =============================================================================
# PHASE 5: Site-Specific Patterns
# =============================================================================

cat("\n========== PHASE 5: Site-Specific Patterns ==========\n\n")

for (sp in c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")) {
  cat("\n", sp, ":\n")
  sp_samples <- rownames(d$sample_table[d$sample_table$SPP == sp, ])
  sp_endo_otu <- endo_otu[, colnames(endo_otu) %in% sp_samples, drop = FALSE]
  sp_endo_otu <- sp_endo_otu[rowSums(sp_endo_otu) > 0, , drop = FALSE]
  
  if (nrow(sp_endo_otu) < 2) {
    cat("  Only", nrow(sp_endo_otu), "ASV(s) - skipping\n"); next
  }
  
  sp_meta <- d$sample_table[sp_samples, ]
  sp_endo_t <- t(sp_endo_otu)
  sp_endo_t <- sp_endo_t[rownames(sp_endo_t) %in% rownames(sp_meta), ]
  sp_meta_matched <- sp_meta[rownames(sp_endo_t), ]
  
  if (length(unique(sp_meta_matched$Site)) > 1) {
    bc_dist <- vegdist(sp_endo_t, method = "bray")
    perm_result <- adonis2(bc_dist ~ Site, data = sp_meta_matched, permutations = 999)
    cat("  PERMANOVA (Endo by Site): R2 =", round(perm_result$R2[1], 4),
        ", p =", perm_result$`Pr(>F)`[1], "\n")
  }
}

# Remove individual annotate tags and adjust tree
p_tree_fig5 <- ggtree(tree_bs, layout = "rectangular", size = 0.5, color = "grey30") %<+%
  tip_data +
  geom_tippoint(aes(color = Host_Group, shape = Source), size = 3, stroke = 0.5) +
  geom_tiplab(aes(label = Short_Label), size = 4.5, offset = 0.005, fontface = "bold") +
  geom_nodelab(aes(label = ifelse(!is.na(as.numeric(label)) & as.numeric(label) >= 70,
                                  round(as.numeric(label)), "")),
               size = 2.5, hjust = -0.2, vjust = -0.5, color = "black") +
  scale_color_manual(
    values = tree_colors, name = "Dominant Host",
    breaks = c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata", "Reference"),
    labels = c(expression(italic("B. asbestinum")), expression(italic("E. caribaeorum")),
               expression(italic("E. flexuosa")), expression(italic("M. muricata")),
               "Reference type strain")
  ) +
  scale_shape_manual(values = c("This study" = 16, "Type strain" = 17), name = "Source") +
  geom_treescale(x = max(fortify(tree_bs)$x) * 0.8, y = -2, 
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

# Haplotype network WITHOUT legend (tree legend covers both)
p_haplo_fig5 <- ggplot() +
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

# Combine with patchwork
fig5 <- p_tree_fig5 + p_haplo_fig5 +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave(file.path(endo_output, "Figure5_Endozoicomonas_combined.png"),
       fig5, width = 14, height = 14, dpi = 300, bg = "white")
ggsave(file.path(endo_output, "Figure5_Endozoicomonas_combined.pdf"),
       fig5, width = 20, height = 14, bg = "white")

cat("Figure 5 saved.\n")
print(fig5)