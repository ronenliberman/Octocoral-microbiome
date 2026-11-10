# =============================================================================
# 00_setup.R — Shared setup for the third-round editor revision (FEMSEC-26-04-0109)
# =============================================================================
# Purpose:  Paths, packages, fixed species palette, the standard dataset
#           filtering block, validation thresholds, and helper functions used
#           by scripts 01–08 in this folder.
# Author:   Ronen Liberman
# Date:     September 2026
#
# Every other script in this folder starts with:
#   source(file.path(script_dir, "00_setup.R"))
# =============================================================================

# ---- PACKAGES ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(microeco)
  library(file2meco)
  library(qiime2R)
  library(tidyverse)
  library(magrittr)
  library(Biostrings)
  library(ape)
  library(ggtree)
  library(writexl)
})

set.seed(123)

# ---- PATHS -------------------------------------------------------------------
# Every path is relative to the repository root, which here::here() finds by
# the .Rproj file. See HELPER.md for the layout.

rev_dir     <- here::here()
script_dir  <- file.path(rev_dir, "script", "Editor_round3")
raw_dir     <- file.path(rev_dir, "raw files")   # QIIME2 outputs (read-only)
blast_dir   <- file.path(rev_dir, "blast")       # saved NCBI BLAST results
old_out_dir <- file.path(rev_dir, "output")      # outputs of the round-2 scripts in script/
supp_dir    <- file.path(rev_dir, "output", "Supplementary")

# New outputs for this round go to their own folder, so nothing from the
# submitted R1 version is overwritten.
out_dir <- file.path(rev_dir, "output", "Editor_round3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tree inputs and results: *_for_tree.fasta and *_reference_names.csv (written
# by 03/04), the MAFFT alignments and the IQ-TREE output (from 11, or the Mac step).
# The name mac_dir is historical: the first alignments and trees were run on a Mac.
mac_dir <- file.path(rev_dir, "phylogeny", "bacteria")
dir.create(mac_dir, showWarnings = FALSE, recursive = TRUE)

# ---- SPECIES PALETTE (single source of truth) --------------------------------
# Matches Fig 2. NOTE: the submitted Fig 3 had M. muricata and E. flexuosa
# swapped; every script in this folder uses this vector instead.

species_colors <- c(
  "B. asbestinum"  = "#E41A1C",
  "E. caribaeorum" = "#FF7F00",
  "E. flexuosa"    = "#984EA3",
  "M. muricata"    = "#4DAF4A"
)
species_order_levels <- c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")

# ---- VALIDATION THRESHOLDS (Yarza et al. 2014, doi:10.1038/nrmicro3330) ------
# Median 16S rRNA identity between members of the same taxon (full-length).
# Used here as approximate cut-offs for ~253 bp V4 reads.

PIDENT_FAMILY <- 86.5   # family  — minimum to accept the family assignment
PIDENT_GENUS  <- 94.5   # genus   — minimum to call an ASV by genus name

# Tree tip selection
MIN_SAMPLES_TREE <- 5   # validated ASVs detected in >= 5 samples (core ASVs always kept)

# Fig 3 row filter
MIN_MEAN_RELABUND_FIG3 <- 1   # % mean relative abundance within the host

# Mycoplasma host assignment: sample counts only when Mycoplasma > 1% of reads
MYCO_MIN_FRACTION <- 0.01

# NCBI family names treated as equivalent to the SILVA family.
# SILVA 138 (the release our classifier was trained on; see qza provenance)
# still uses the pre-2018 Mollicutes taxonomy; NCBI has split
# Mycoplasmataceae (Gupta et al. 2018). Check these before running 01.
family_equivalents <- list(
  Endozoicomonadaceae = c("Endozoicomonadaceae"),
  Mycoplasmataceae    = c("Mycoplasmataceae", "Metamycoplasmataceae",
                          "Mycoplasmoidaceae")
)

# Core ASVs (hashes) — fixed labels kept from the R1 submission so the text
# and Table S8 stay valid.
core_endo_map <- c(
  "238ab607aedd1bd3b1b75b03f61331c3" = "Endo_1",
  "2f40b4eaa74d3e96797cbb513a5f43bd" = "Endo_2",
  "47ce1658b2ee5601bd633f5f2d014adf" = "Endo_3",
  "5656aae49e1f57ea6217094b1969e8aa" = "Endo_4",
  "69d03d98367f96a5ec36f2c90e977fad" = "Endo_5",
  "a80ddcda6c71332a7a6269bdce4dee87" = "Endo_6"
)
# Myco_1–Myco_10 are the 10 most abundant Mycoplasma ASVs by total reads
# (as in Mycoplasma_Phylogenetic_Analysis_UPDATED.R); taken from the R1
# master table in 02_stable_ASV_IDs.R.

# =============================================================================
# REFERENCE SEQUENCES FOR THE TREES — NCBI RefSeq type-strain 16S (NR_)
# =============================================================================
# The R1 trees used GenBank accessions, several of which were NOT the labelled
# species (checked 1 Oct 2026 with NCBI esummary):
#   MG456764 "E. lisbonensis" = Vibrio scophthalmi
#   KY770540 "E. acroporae"   = Roseobacter sp.;  KY770538 "E. ascidiicola" = Roseobacter sp.
#   NR_025989 "M. iowae" = Mesomycoplasma ovipneumoniae; NR_113669 "M. todarodis" = Sinorhizobium fredii
#   NR_113899 "M. insons" = Enterococcus avium; NR_025949 "Acholeplasma laidlawii" = Bordetella bronchiseptica
# All references are now RefSeq NR_ records, looked up by organism name.

endo_references <- tribble(
  ~Accession,    ~Label,                         ~Outgroup,
  "NR_109684.2", "E. euniceicola EF212T",         FALSE,
  "NR_109685.1", "E. gorgoniicola PS125T",        FALSE,
  "NR_041264.1", "E. elysicola MKT110T",          FALSE,
  "NR_114318.1", "E. numazuensis HC50T",          FALSE,
  "NR_116609.1", "E. montiporae CL-33T",          FALSE,
  "NR_134024.1", "E. atrinae WP70T",              FALSE,
  "NR_158127.1", "E. acroporae Acr-14T",          FALSE,
  "NR_146693.1", "E. ascidiicola AVMART05T",      FALSE,
  "NR_199442.1", "E. lisbonensis NE40T",          FALSE,
  "NR_116386.1", "Kistimonas asteriae KMD 001T",  TRUE
)

myco_references <- tribble(
  ~Accession,    ~Label,                           ~Outgroup,
  "NR_179594.1", "M. todarodis 5HT",               FALSE,
  "NR_044669.2", "Malacoplasma iowae 695T",        FALSE,
  "NR_115812.1", "M. insons I17P1T",               FALSE,
  "NR_074448.2", "Acholeplasma laidlawii PG-8AT",  TRUE
)

# 515F GTGYCAGCMGCCGCGGTAA ... 806R reverse complement ATTAGAWACCCBDGTHGTCC.
# Keep the sequence between the primers, the same region as the ASVs (~253 bp).
trim_to_v4_amplicon <- function(full_seq) {
  f <- regexpr("GTG[CT]CAGC[AC]GCCGCGGTAA", full_seq, ignore.case = TRUE)
  if (f == -1) f <- regexpr("CAGC[AC]GCCGCGGTAA", full_seq, ignore.case = TRUE)
  if (f == -1) return(NA_character_)
  start <- f + attr(f, "match.length")
  rest  <- substr(full_seq, start, nchar(full_seq))
  r <- regexpr("ATTAGA[AT]ACCC", rest, ignore.case = TRUE)
  if (r == -1 || r < 200 || r > 300) return(NA_character_)
  toupper(substr(rest, 1, r - 1))
}

# Download (efetch, one batch), check the organism in the record title matches
# the label's genus, trim to the V4 amplicon, and return a named DNAStringSet
# (IQ-TREE-safe names "REF_...") plus a table of what was used.
get_reference_v4 <- function(ref_tbl) {
  txt <- httr::content(httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                 query = list(db = "nuccore", id = paste(ref_tbl$Accession, collapse = ","),
                                              rettype = "fasta", retmode = "text")),
                       "text", encoding = "UTF-8")
  recs <- strsplit(txt, "\n(?=>)", perl = TRUE)[[1]]
  recs <- recs[nzchar(trimws(recs))]
  parsed <- tibble(
    header   = sub("\n.*", "", recs),
    Accession = sub("^>(\\S+).*", "\\1", sub("\n.*", "", recs)),
    Sequence = gsub("\\s", "", sub("^[^\n]*\n", "", recs))
  )
  out <- ref_tbl %>%
    left_join(parsed, by = "Accession") %>%
    mutate(
      Genus_in_label  = case_when(grepl("^E\\. ", Label) ~ "Endozoicomonas",
                                  grepl("^M\\. ", Label) ~ "Mycoplasma",
                                  TRUE ~ word(Label, 1)),
      Record_title    = sub("^>\\S+\\s+", "", header),
      V4              = map_chr(Sequence, trim_to_v4_amplicon),
      V4_length       = nchar(V4),
      Safe            = paste0("REF_", gsub("[^A-Za-z0-9]+", "_", Label))
    )
  # the label's genus must appear in the NCBI record title
  out$Title_check <- mapply(function(g, t) grepl(g, t, ignore.case = TRUE),
                            out$Genus_in_label, out$Record_title)
  if (any(!out$Title_check)) {
    print(out %>% filter(!Title_check) %>% dplyr::select(Accession, Label, Record_title))
    stop("Reference title does not match its label — check the accession.")
  }
  if (any(is.na(out$V4))) {
    print(out %>% filter(is.na(V4)) %>% dplyr::select(Accession, Label))
    stop("Could not find the 515F/806R sites in these references.")
  }
  seqs <- DNAStringSet(out$V4); names(seqs) <- out$Safe
  list(seqs = seqs, table = out %>% dplyr::select(Accession, Label, Outgroup, Record_title, V4_length, Safe))
}

# =============================================================================
# DATASET — standard filtering block (identical to the other revision scripts)
# =============================================================================

load_octocoral_dataset <- function(include_seawater = FALSE) {
  cat("Loading QIIME2 artifacts (PEE2)...\n")
  data <- qiime2meco(
    file.path(raw_dir, "table-NoEMC-PEE2.qza"),
    taxonomy_table = file.path(raw_dir, "taxonomy-PEE2.qza"),
    sample_table   = file.path(raw_dir, "metadata_all_edited.csv"),
    auto_tidy      = TRUE
  )
  data$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

  keep_spp <- c("BAST", "ECAR", "MMUR", "EFLE", if (include_seawater) "SW")
  d <- clone(data)
  d$sample_table %<>% subset(SPP %in% keep_spp)
  d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
  d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

  samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024", "EFLE_DC3_3_2024")
  d$sample_table <- d$sample_table %>% filter(!rownames(.) %in% samples_to_remove)

  species.names <- list(
    "SW"   = "Seawater",
    "BAST" = "B. asbestinum",
    "ECAR" = "E. caribaeorum",
    "MMUR" = "M. muricata",
    "EFLE" = "E. flexuosa"
  )
  d$sample_table <- d$sample_table %>% mutate(SPP = recode(SPP, !!!species.names))
  d$tidy_dataset()
  cat("Filtered dataset:", nrow(d$sample_table), "samples,",
      nrow(d$otu_table), "ASVs\n")
  d
}

load_rep_seqs <- function() {
  readDNAStringSet(file.path(raw_dir, "dna-sequences-PEE2.fasta"))
}

# SILVA classifier confidence per ASV (Confidence column of taxonomy-PEE2.qza)
load_silva_confidence <- function() {
  tx <- read_qza(file.path(raw_dir, "taxonomy-PEE2.qza"))$data
  tibble(ASV_ID = tx$Feature.ID, SILVA_confidence = as.numeric(tx$Confidence))
}

# =============================================================================
# HELPERS
# =============================================================================

strip_prefix <- function(x) sub("^[a-z]__", "", x)

# ---- Per-ASV summary (samples, reads, hosts, dominant host) ------------------
# Dominant host = host species with the most samples containing the ASV.
# If min_group_fraction is set, a sample only counts when the summed reads of
# `group_asvs` (e.g. all Mycoplasma ASVs) exceed that fraction of the sample.

build_asv_master <- function(d, asv_ids, group_asvs = NULL, min_group_fraction = NULL) {
  otu <- d$otu_table
  asv_ids <- intersect(asv_ids, rownames(otu))
  sample_totals <- colSums(otu)
  meta <- d$sample_table

  keep_sample <- rep(TRUE, ncol(otu)); names(keep_sample) <- colnames(otu)
  if (!is.null(min_group_fraction)) {
    grp <- colSums(otu[intersect(group_asvs, rownames(otu)), , drop = FALSE])
    keep_sample <- (grp / sample_totals) > min_group_fraction
  }

  bind_rows(lapply(asv_ids, function(asv) {
    counts  <- unlist(otu[asv, ])
    present <- names(counts)[counts > 0]
    if (length(present) == 0) return(NULL)
    present_kept <- present[keep_sample[present]]

    hosts_all <- table(meta[present, "SPP"])
    hosts_kept <- table(meta[present_kept, "SPP"])
    dom <- if (length(hosts_kept) == 0 || max(hosts_kept) == 0) "Low abundance" else
      names(hosts_kept)[which.max(hosts_kept)]

    tibble(
      ASV_ID            = asv,
      N_Samples         = length(present),
      Total_Reads       = sum(counts),
      Mean_RelAbund_Pct = round(mean((counts / sample_totals)[present]) * 100, 4),
      Host_Species      = paste(sort(names(hosts_all)), collapse = "; "),
      N_Host_Species    = length(hosts_all),
      N_Sites           = n_distinct(meta[present, "Site"]),
      Dominant_Host     = dom
    )
  })) %>% arrange(desc(Total_Reads))
}

# ---- Alignment statistics ----------------------------------------------------
# total columns, columns with no gap/ambiguity in ANY sequence ("unambiguously
# aligned"), and variable / parsimony-informative columns among those.

alignment_stats <- function(aln_path, label) {
  aln <- as.character(read.dna(aln_path, format = "fasta"))
  aln <- toupper(aln)
  n_seq <- nrow(aln); n_col <- ncol(aln)
  gap_free <- apply(aln, 2, function(col) all(col %in% c("A", "C", "G", "T")))
  sub <- aln[, gap_free, drop = FALSE]
  variable <- apply(sub, 2, function(col) length(unique(col)) > 1)
  pars_inf <- apply(sub, 2, function(col) sum(table(col) >= 2) >= 2)
  tibble(
    Alignment                      = label,
    N_sequences                    = n_seq,
    Alignment_length               = n_col,
    Unambiguously_aligned_positions = sum(gap_free),
    Variable_positions             = sum(variable),
    Parsimony_informative_positions = sum(pars_inf)
  )
}

# Pull the site summary and model from an IQ-TREE .iqtree report
iqtree_summary <- function(iqtree_file) {
  txt <- readLines(iqtree_file, warn = FALSE)
  grab <- function(pattern) {
    hit <- grep(pattern, txt, value = TRUE)
    if (length(hit) == 0) NA_character_ else trimws(hit[1])
  }
  tibble(
    File            = basename(iqtree_file),
    Input_line      = grab("^Input data:"),
    Parsimony_line  = grab("parsimony-informative"),
    Best_model      = sub(".*:\\s*", "", grab("^Best-fit model according to BIC"))
  )
}

# ---- ML tree plot (tree only, every tip labelled) ------------------------------
# tip_data columns required: label, Tip_Label (display text), Is_Ref (logical),
# Host_Group (species name, "Low abundance" or "Reference"), Is_Core (logical)

# Branch support: IQ-TREE node labels are "SH-aLRT/UFBoot". A node is supported
# when SH-aLRT >= 80% AND UFBoot >= 95% (Guindon et al. 2010; Hoang et al. 2018);
# nodes meeting only one of the two criteria are drawn grey, others not drawn.
plot_ml_tree <- function(tree, tip_data, outgroup_pattern,
                         point_size = 6, label_size = 5, scale_width = 0.01,
                         extra_colors = c(), sh_cut = 80, uf_cut = 95,
                         legend_y = 0.5, layout = list()) {
  # Layout options (v2 figures, 19 Oct 2026); the defaults reproduce the _FIXED figures.
  #   xlim_mult: x range as a multiple of the deepest tip; legend_x/legend_just:
  #   legend anchor (npc); legend_scale: enlarges legend text, keys and symbols;
  #   outgroup_above: outgroup name written above its own branch, ending at the tip;
  #   legend_bottom: legends in rows under the tree (host_nrow / source_nrow rows,
  #   titles on the left, key_gap_cm between entries; see plot_ml_tree_support_values);
  #   legend_box (v3, 20 Oct 2026): the same rows as an outlined inset in the tree's
  #   white space, anchored at (legend_x, legend_y) with legend_just
  lo <- modifyList(list(xlim_mult = 1.45, legend_x = 0.98, legend_just = c(1, legend_y),
                        legend_scale = 1, outgroup_above = FALSE, legend_bottom = FALSE,
                        legend_box = FALSE, host_nrow = NULL, source_nrow = NULL,
                        key_gap_cm = 0.5), layout)
  rows <- lo$legend_bottom || lo$legend_box

  options(ignore.negative.edge = TRUE)
  out_tip <- grep(outgroup_pattern, tree$tip.label, value = TRUE)
  # edgelabel = TRUE: SH-aLRT/UFBoot belong to splits and must move with their branches when re-rooting (fixed 14 Oct 2026, audit A9)
  tree <- root(tree, outgroup = out_tip, resolve.root = TRUE, edgelabel = TRUE)
  tree <- ladderize(tree)

  both_lab <- paste0("SH-aLRT ≥ ", sh_cut, "% and UFBoot ≥ ", uf_cut, "%")
  one_lab  <- "One criterion met"
  node_df <- fortify(tree) %>%
    filter(!isTip) %>%
    mutate(sh = suppressWarnings(as.numeric(sub("/.*", "", label))),
           uf = suppressWarnings(as.numeric(sub(".*/", "", label))),
           sh_ok = !is.na(sh) & sh >= sh_cut,
           uf_ok = !is.na(uf) & uf >= uf_cut) %>%
    filter(grepl("/", label), sh_ok | uf_ok) %>%
    mutate(bs_tier = factor(ifelse(sh_ok & uf_ok, both_lab, one_lab),
                            levels = c(one_lab, both_lab)))
  cat("Branch support symbols:", sum(node_df$bs_tier == both_lab), "nodes with both criteria,",
      sum(node_df$bs_tier == one_lab), "with one criterion\n")

  # Italic for type strains, bold for core ASVs, plain for the rest
  tip_data <- tip_data %>%
    mutate(Tip_Label_parsed = case_when(
      Is_Ref  ~ paste0("italic('", Tip_Label, "')"),
      Is_Core ~ paste0("bold('", Tip_Label, "')"),
      TRUE    ~ paste0("'", Tip_Label, "'")
    ))

  tree_colors <- c(species_colors, "Reference" = "black",
                   "Low abundance" = "grey70", extra_colors)
  # "Reference" is shown by the Source (triangle) legend, not repeated under Host
  present_groups <- setdiff(intersect(names(tree_colors), unique(tip_data$Host_Group)), "Reference")
  legend_strings <- sapply(present_groups, function(g) {
    if (g %in% names(species_colors)) paste0("italic('", g, "')") else
      if (g == "Reference") "'Type strain'" else
        if (g == "Low abundance") "'Trace only (<1%)'" else paste0("'", g, "'")
  })
  legend_labels <- parse(text = legend_strings)

  max_x <- max(fortify(tree)$x)
  ls_ <- lo$legend_scale
  ts_ <- if (is.null(lo$legend_text_scale)) ls_ else lo$legend_text_scale   # legend text/titles only
  tip_data$Label_Above <- lo$outgroup_above & tip_data$label %in% out_tip

  # Legend in rows: pad the narrower title and first-column labels so the Host and
  # Source keys line up in columns (widths measured in pt at print size)
  host_th <- src_th <- theme()
  if (rows) {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
    w_pt <- function(s, face) max(sapply(s, function(z) grid::convertWidth(grid::grobWidth(grid::textGrob(
      z, gp = grid::gpar(fontsize = (if (face == 2) 15 else 14) * ts_, fontface = face))), "pt", valueOnly = TRUE)))
    d_title <- w_pt("Source", 2) - w_pt("Host", 2)
    host_col1 <- present_groups[seq_len(min(length(present_groups), max(1, lo$host_nrow)))]   # filled by column
    d_lab   <- w_pt(host_col1, 3) - w_pt("This study", 1)
    gap_pt  <- lo$key_gap_cm / 2.54 * 72
    host_th <- theme(legend.title = element_text(margin = margin(r = 8 + max(0, d_title))),
                     legend.key.spacing.x = unit(gap_pt + max(0, -d_lab), "pt"))
    src_th  <- theme(legend.title = element_text(margin = margin(r = 8 + max(0, -d_title))),
                     legend.key.spacing.x = unit(gap_pt + max(0, d_lab), "pt"))
  }
  rows_theme <- theme(
    legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "left",
    legend.title.position = "left", legend.title = element_text(size = 15 * ts_, face = "bold", hjust = 0),
    legend.key.spacing.x = unit(lo$key_gap_cm, "cm")
  )
  legend_theme <- if (lo$legend_box) rows_theme + theme(
    legend.position = "inside", legend.position.inside = c(lo$legend_x, legend_y),
    legend.justification.inside = lo$legend_just,
    legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    legend.box.margin = margin(6, 12, 6, 12),
    legend.key.height = unit(0.6 * ls_, "cm"),   # compact rows
    legend.spacing.y = unit(4, "pt"),            # Host block to Source block
    legend.margin = margin(2, 2, 2, 2)
  ) else if (lo$legend_bottom) rows_theme + theme(
    legend.position = "bottom"
  ) else theme(
    legend.position  = c(lo$legend_x, legend_y),
    legend.justification = lo$legend_just
  )

  ggtree(tree, size = 0.7, color = "grey25") %<+% tip_data +
    geom_tippoint(aes(color = Host_Group, shape = Is_Ref), size = point_size, stroke = 0.6) +
    geom_tiplab(aes(label = Tip_Label_parsed, subset = !Label_Above), parse = TRUE,
                size = label_size, offset = max_x * 0.02, align = FALSE) +
    geom_text2(aes(x = x - max_x * 0.03, y = y + 0.12, label = Tip_Label_parsed,
                  subset = isTip & Label_Above),
              parse = TRUE, size = label_size, hjust = 1, vjust = 0) +
    geom_point2(data = node_df, aes(x = x, y = y, fill = bs_tier),
                shape = 22, size = 3.5, color = "black", stroke = 0.4,
                inherit.aes = FALSE) +
    scale_fill_manual(values = setNames(c("grey60", "black"), c(one_lab, both_lab)),
                      name = "Branch support", drop = FALSE) +
    scale_color_manual(values = tree_colors, breaks = present_groups,
                       labels = legend_labels, name = "Host") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                       labels = c("FALSE" = "This study", "TRUE" = "Type strain"),
                       name = "Source") +
    # scale bar in its own row below the tree, so it never collides with labels
    geom_treescale(x = 0, y = -0.6, width = scale_width,
                   fontsize = 5, linesize = 1, offset = 0.5) +
    xlim(0, max_x * lo$xlim_mult) +     # room for labels, tree uses the page width
    ylim(-1.2, Ntip(tree) + 0.5) +
    theme_tree() +
    theme(
      legend.title     = element_text(size = 15 * ts_, face = "bold"),
      legend.text      = element_text(size = 14 * ts_),
      legend.key.size  = unit(0.8 * ls_, "cm"),
      legend.background = element_blank(),
      plot.margin      = margin(10, 10, 10, 10)
    ) +
    legend_theme +
    guides(color = guide_legend(override.aes = list(size = 6 * ls_), nrow = lo$host_nrow, order = 1, theme = host_th),
           shape = guide_legend(override.aes = list(size = 6 * ls_), nrow = lo$source_nrow, order = 2, theme = src_th),
           fill  = guide_legend(override.aes = list(size = 5), order = 3))
}

# ---- ML tree with "SH-aLRT/UFBoot" values at the nodes (final Figs 5-6) --------
# Adopted 4 Oct 2026 after the visual trials in 13_TRIAL_tree_UFBoot_labels.R
# (see Editor_Revision_Fix_Plan.md, "Branch-support display"). Builds on
# plot_ml_tree() (squares switched off) and prints SH-aLRT/UFBoot for every node
# with UFBoot >= uf_min. SH-aLRT is rounded but never shown as 100 when < 100.
# Positions are chosen by place_support_values(), a collision checker in inches
# at print size: left of the node above/below its branch (convention), or inside
# the node's own fork; a value must touch no branch, tip symbol, tip label or
# other value, stay inside the plot and be clearly nearer its own node (<= 75%
# of the distance to any other node). A stacked "SH/" over "UF" form is used
# only if no single-line position is free. fig_w/fig_h must equal the ggsave size.
place_support_values <- function(p, nodes, fig_w, fig_h, value_size, tip_label_size = 5,
                                 pad_in = 0.025, margin_in = 10 / 72, xlim_mult = 1.45) {
  d     <- p$data
  max_x <- max(d$x)
  xr    <- c(0, max_x * xlim_mult)                 # xlim() in plot_ml_tree()
  yr    <- c(-1.2, sum(d$isTip) + 0.5)             # ylim() in plot_ml_tree()
  sx <- (fig_w - 2 * margin_in) / diff(xr)          # inches per x unit
  sy <- (fig_h - 2 * margin_in) / diff(yr)          # inches per y unit
  X <- function(x) (x - xr[1]) * sx; Y <- function(y) (y - yr[1]) * sy

  grDevices::pdf(NULL, width = fig_w, height = fig_h); on.exit(grDevices::dev.off())
  txt_w <- function(s, size) sapply(s, function(z) grid::convertWidth(
    grid::grobWidth(grid::textGrob(z, gp = grid::gpar(fontsize = size * .pt))), "inches", valueOnly = TRUE))
  h_val <- value_size * .pt / 72                    # text height (in)

  # obstacles (inches): branch segments, tip symbols, tip labels
  par <- d[match(d$parent, d$node), ]
  segs <- bind_rows(
    tibble(x1 = X(par$x), x2 = X(d$x), y1 = Y(d$y), y2 = Y(d$y)),                 # horizontal
    d %>% filter(!isTip) %>% rowwise() %>%                                         # vertical
      mutate(lo = min(d$y[d$parent == node & d$node != node]), hi = max(d$y[d$parent == node & d$node != node])) %>%
      ungroup() %>% transmute(x1 = X(x), x2 = X(x), y1 = Y(lo), y2 = Y(hi))) %>%
    filter(is.finite(y1), is.finite(y2))
  r_pt  <- 0.10                                                                   # tip symbol radius (in)
  tips  <- d %>% filter(isTip)
  boxes <- bind_rows(
    tibble(xmin = X(tips$x) - r_pt, xmax = X(tips$x) + r_pt, ymin = Y(tips$y) - r_pt, ymax = Y(tips$y) + r_pt),
    tibble(xmin = X(tips$x + max_x * 0.02), ymin = Y(tips$y) - tip_label_size * .pt / 144,
           ymax = Y(tips$y) + tip_label_size * .pt / 144) %>%
      mutate(xmax = xmin + txt_w(as.character(tips$label), tip_label_size) * 1.15))

  seg_hits <- function(b) sum(pmax(0, pmin(segs$x2, b["xmax"]) - pmax(segs$x1, b["xmin"]) + (segs$x1 == segs$x2) * 1e-9) > 0 &
                              pmax(0, pmin(pmax(segs$y1, segs$y2), b["ymax"]) - pmax(pmin(segs$y1, segs$y2), b["ymin"]) + (segs$y1 == segs$y2) * 1e-9) > 0)
  box_hits <- function(b, B) if (!nrow(B)) 0 else sum(B$xmin < b["xmax"] & B$xmax > b["xmin"] & B$ymin < b["ymax"] & B$ymax > b["ymin"])
  inner <- d %>% filter(!isTip) %>% transmute(nx2 = X(x), ny2 = Y(y), node)

  gap_x <- 0.06; gap_y <- 0.03
  placed <- tibble(xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric())
  nodes$lab_x <- NA_real_; nodes$lab_y <- NA_real_; nodes$where <- NA_character_
  for (i in order(-nodes$y)) {
    nx <- X(nodes$x[i]); ny <- Y(nodes$y[i])
    parts <- strsplit(nodes$txt[i], "/")[[1]]
    forms <- list(inline = list(txt = nodes$txt[i], w = unname(txt_w(nodes$txt[i], value_size)), h = h_val, pen = 0),
                  stack  = list(txt = paste0(parts[1], "/\n", parts[2]),
                                w = max(unname(txt_w(c(paste0(parts[1], "/"), parts[2]), value_size))),
                                h = 2 * h_val * 0.95, pen = 0.5))
    cands <- list()
    kids <- d %>% filter(parent == nodes$node[i], node != nodes$node[i])           # (a) inside own fork
    if (nrow(kids) >= 2) {
      lo <- Y(min(kids$y)); hi <- Y(max(kids$y))
      for (fm in names(forms)) {
        f <- forms[[fm]]; w <- f$w + 2 * pad_in; hh <- f$h
        y_from <- lo + hh / 2 + pad_in + 0.03; y_to <- hi - hh / 2 - pad_in - 0.03
        if (y_from > y_to) next
        for (yc in seq(y_from, y_to, by = 0.02)) {
          xr_ <- nx + 0.05 + w
          cands[[length(cands) + 1]] <- list(b = c(xmin = xr_ - w, xmax = xr_, ymin = yc - hh / 2 - pad_in, ymax = yc + hh / 2 + pad_in),
                                             xr = xr_, yc = yc, txt = f$txt, pen = f$pen + 0.3 + abs(yc - ny) * 0.2,
                                             where = paste0(fm, ", inside fork"), fork = TRUE)
        }
      }
    }
    for (fm in names(forms)) for (gx in c(gap_x, 0.025)) for (step in 0:8)        # (b) left of node
      for (side in c("above", "below")) for (lift in c(0, 1)) {
        if (gx != gap_x && step > 0) next
        f   <- forms[[fm]]; w <- f$w + 2 * pad_in; hh <- f$h
        xr_ <- nx - gx - step * 0.06
        dyl <- lift * (h_val * 0.6)
        yc  <- if (side == "above") ny + gap_y + hh / 2 + pad_in + dyl else ny - gap_y - hh / 2 - pad_in - dyl
        cands[[length(cands) + 1]] <- list(b = c(xmin = xr_ - w, xmax = xr_, ymin = yc - hh / 2 - pad_in, ymax = yc + hh / 2 + pad_in),
                                           xr = xr_, yc = yc, txt = f$txt, pen = f$pen + (gx != gap_x) * 0.05,
                                           where = paste0(fm, ", ", side, ifelse(step, paste0(", left ", step * 0.06, " in"), ""),
                                                          ifelse(lift, ", +lift", ""), ifelse(gx != gap_x, ", tight", "")),
                                           fork = FALSE)
      }
    anchor_x <- function(cc) if (isTRUE(cc$fork)) unname(cc$b["xmin"]) else cc$xr   # edge facing the node
    own_closer <- function(cc) {
      dd <- sqrt((inner$nx2 - anchor_x(cc))^2 + (inner$ny2 - cc$yc)^2)
      dd[inner$node == nodes$node[i]] <= 0.75 * min(dd[inner$node != nodes$node[i]])
    }
    score <- unname(sapply(cands, function(cc) {
      clash <- seg_hits(cc$b) + box_hits(cc$b, boxes) + box_hits(cc$b, placed) +
               (cc$b["xmin"] < X(0) + 0.02) * 10 + (!own_closer(cc)) * 5
      clash * 100 + cc$pen + sqrt((anchor_x(cc) - nx)^2 + (cc$yc - ny)^2)
    }))
    k  <- which.min(score); cc <- cands[[k]]
    placed <- bind_rows(placed, as_tibble(as.list(cc$b)))
    nodes$lab_x[i] <- cc$xr / sx + xr[1] - pad_in / sx
    nodes$lab_y[i] <- cc$yc / sy + yr[1]
    nodes$txt[i]   <- cc$txt
    nodes$where[i] <- paste0(cc$where, ifelse(score[k] >= 100, paste0("  ** conflict score ", floor(score[k] / 100)), ""))
  }
  nodes
}

plot_ml_tree_support_values <- function(tree, tip_data, outgroup_pattern, fig_w, fig_h,
                                        legend_y = 0.5, uf_min = 70, value_size = 3.6, tag = "",
                                        layout = list()) {
  p <- plot_ml_tree(tree, tip_data, outgroup_pattern, point_size = 6, label_size = 5,
                    scale_width = 0.01, sh_cut = Inf, uf_cut = Inf, legend_y = legend_y,
                    layout = layout)
  # Bottom legend: the tree is laid out at fig_h - legend_h_in and the legend strip
  # is added below it, so the inch-based placement below matches the printed tree
  leg_h  <- if (isTRUE(layout$legend_bottom)) layout$legend_h_in else 0
  tree_h <- fig_h - leg_h
  nodes <- p$data %>%
    filter(!isTip, grepl("/", label)) %>%
    mutate(sh = suppressWarnings(as.numeric(sub("/.*", "", label))),
           uf = suppressWarnings(as.numeric(sub(".*/", "", label)))) %>%
    filter(!is.na(uf), uf >= uf_min) %>%
    transmute(node, x, y, txt = paste0(ifelse(sh < 100, pmin(round(sh), 99), 100), "/", round(uf)))
  nodes <- place_support_values(p, nodes, fig_w, tree_h, value_size,
                                xlim_mult = if (is.null(layout$xlim_mult)) 1.45 else layout$xlim_mult)
  cat("\n", tag, "SH-aLRT/UFBoot values drawn:", nrow(nodes),
      "| conflicts:", sum(grepl("conflict", nodes$where)), "\n")
  print(as.data.frame(nodes %>% mutate(txt = gsub("\n", " | ", txt)) %>% dplyr::select(node, txt, where)), row.names = FALSE)

  max_x <- max(p$data$x); hx <- max_x * 0.0012; hy <- 0.035                       # white halo
  halo <- tidyr::crossing(nodes, tibble(ox = c(-1, 0, 1, -1, 1, -1, 0, 1), oy = c(-1, -1, -1, 0, 0, 1, 1, 1))) %>%
    mutate(hx_ = lab_x + ox * hx, hy_ = lab_y + oy * hy)
  p <- p +
    geom_text(data = halo, aes(x = hx_, y = hy_, label = txt), inherit.aes = FALSE,
              size = value_size, colour = "white", hjust = 1, vjust = 0.5, lineheight = 0.85) +
    geom_text(data = nodes, aes(x = lab_x, y = lab_y, label = txt), inherit.aes = FALSE,
              size = value_size, colour = "black", hjust = 1, vjust = 0.5, lineheight = 0.85) +
    guides(fill = "none")
  if (leg_h == 0) return(p)
  leg <- cowplot::get_plot_component(p, "guide-box-bottom", return_all = FALSE)
  cowplot::plot_grid(p + theme(legend.position = "none"), leg, ncol = 1,
                     rel_heights = c(tree_h, leg_h))
}

# ---- Saved tree + tip table, for redrawing Figs 5-6 without re-running 03/04 -----
# Reads the IQ-TREE tree (mac_dir) and the tip tables 03/04 wrote; labels and
# host groups follow the same rules as 03/04 (dagger = 04's † for unverified ASVs).
tree_tip_data <- function(prefix, status_csv, dagger = FALSE) {
  tree    <- read.tree(file.path(mac_dir, paste0(prefix, "_iqtree.treefile")))
  tips    <- read.csv(file.path(out_dir, status_csv), stringsAsFactors = FALSE)
  ref_map <- read.csv(file.path(mac_dir, paste0(prefix, "_reference_names.csv")), stringsAsFactors = FALSE)
  td <- tibble(label = tree$tip.label) %>%
    mutate(Is_Ref = grepl("^REF_", label)) %>%
    left_join(tips %>% dplyr::select(Stable_ID, Dominant_Host, Is_Core, Verified),
              by = c("label" = "Stable_ID")) %>%
    left_join(ref_map %>% dplyr::select(Display, Safe), by = c("label" = "Safe")) %>%
    mutate(Tip_Label  = case_when(Is_Ref ~ Display,
                                  dagger & !Verified ~ paste0(label, "†"),
                                  TRUE ~ label),
           Host_Group = ifelse(Is_Ref, "Reference", Dominant_Host),
           Is_Core    = coalesce(Is_Core, FALSE))
  stopifnot(!any(is.na(td$Host_Group)))
  list(tree = tree, tip_data = td)
}

cat("00_setup.R loaded. Outputs go to:", out_dir, "\n")
