# =============================================================================
# COMBINED 4-SPECIES COMMUNITY COMPOSITION — UNIFIED LEGEND
# Dominant genera defined as >20% in any sample (excluding uncl. Bacteria)
# Faceted by Species and Year
# =============================================================================

library(tidyverse)
library(microeco)
library(RColorBrewer)
library(cowplot)

output_dir <- file.path("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Output")


# =============================================================================
# ASSESS 100% Brevifollis SAMPLE IN E. flexuosa
# =============================================================================

cat("\n=== Assessing 100% Brevifollis sample ===\n")

# Find the sample(s) where Brevifollis = 100%
efle_data <- clone(mt_TSS)
efle_data$sample_table <- subset(efle_data$sample_table, SPP == "E. flexuosa")
efle_data$tidy_dataset()

otu <- efle_data$otu_table
tax <- efle_data$tax_table

# Get Brevifollis ASVs
brev_asvs <- rownames(tax)[grepl("Brevifollis", tax$Genus)]
cat("Brevifollis ASVs:", length(brev_asvs), "\n")

# Calculate relative abundance of Brevifollis per sample
brev_rel <- colSums(otu[brev_asvs, , drop = FALSE]) / colSums(otu) * 100

# Find samples where Brevifollis > 90%
suspect_samples <- names(brev_rel[brev_rel > 90])
cat("Samples with >90% Brevifollis:", length(suspect_samples), "\n")
cat("Sample IDs:", paste(suspect_samples, collapse = ", "), "\n\n")

# Now check read counts in the UNRAREFIED data
for (samp in suspect_samples) {
  
  # From unrarefied data
  if (samp %in% colnames(d$otu_table)) {
    total_reads_unrarefied <- sum(d$otu_table[, samp])
  } else {
    total_reads_unrarefied <- NA
  }
  
  # From rarefied data
  if (samp %in% colnames(mt_rarefied$otu_table)) {
    total_reads_rarefied <- sum(mt_rarefied$otu_table[, samp])
  } else {
    total_reads_rarefied <- NA
  }
  
  # TSS total (should be ~1)
  total_tss <- sum(otu[, samp])
  brev_pct <- brev_rel[samp]
  n_asvs_detected <- sum(otu[, samp] > 0)
  
  cat("Sample:", samp, "\n")
  cat("  Brevifollis rel. abund:", round(brev_pct, 2), "%\n")
  cat("  Total reads (unrarefied):", total_reads_unrarefied, "\n")
  cat("  Total reads (rarefied):", total_reads_rarefied, "\n")
  cat("  ASVs detected (TSS):", n_asvs_detected, "\n")
  cat("  Metadata:\n")
  print(efle_data$sample_table[samp, ])
  cat("\n")
}

# Compare to other E. flexuosa samples
cat("--- Comparison to other E. flexuosa samples ---\n")
efle_unrar <- clone(d)
efle_unrar$sample_table <- subset(efle_unrar$sample_table, SPP == "E. flexuosa")
efle_unrar$tidy_dataset()

all_reads <- colSums(efle_unrar$otu_table)
all_asvs <- colSums(efle_unrar$otu_table > 0)

cat("E. flexuosa read count summary (unrarefied):\n")
cat("  Mean:", round(mean(all_reads)), "\n")
cat("  Median:", round(median(all_reads)), "\n")
cat("  Min:", min(all_reads), "| Max:", max(all_reads), "\n")
cat("  Suspect sample reads:", all_reads[suspect_samples], "\n\n")

cat("E. flexuosa ASV count summary (unrarefied):\n")
cat("  Mean:", round(mean(all_asvs)), "\n")
cat("  Median:", round(median(all_asvs)), "\n")
cat("  Min:", min(all_asvs), "| Max:", max(all_asvs), "\n")
cat("  Suspect sample ASVs:", all_asvs[suspect_samples], "\n")


# =============================================================================
# DOMINANT GENERA LIST (from >20% analysis, excluding uncl. Bacteria)
# =============================================================================

# All genera that exceeded 20% relative abundance in at least one sample
# across any of the four species. "uncl. Bacteria" excluded.
top_genera_all <- c(
  "Endozoicomonas", "Gammaproteobacteria", "Rhodobacteraceae", "Oleiphilus",
  "Alphaproteobacteria", "Rhizobiaceae", "Mycobacterium",
  "Mycoplasma", "Bacilli", "Tistlia", "Defluviicoccus",
  "Methyloligellaceae", "Labrenzia", "Spiroplasma", "Blastopirellula",
  "SAR116_clade", "Candidatus_Hepatoplasma", "Cyclobacteriaceae",
  "Spirochaeta_2"
)

# Display names with rank prefixes for non-genus groups
top_genera_display <- c(
  "Endozoicomonas"          = "Endozoicomonas",
  "Gammaproteobacteria"     = "c_Gammaproteobacteria",
  "Rhodobacteraceae"        = "f_Rhodobacteraceae",
  "Oleiphilus"              = "Oleiphilus",
  "Alphaproteobacteria"     = "c_Alphaproteobacteria",
  "Rhizobiaceae"            = "f_Rhizobiaceae",
  "Mycobacterium"           = "Mycobacterium",
  "Mycoplasma"              = "Mycoplasma",
  "Bacilli"                 = "c_Bacilli",
  "Tistlia"                 = "Tistlia",
  "Defluviicoccus"          = "Defluviicoccus",
  "Methyloligellaceae"      = "f_Methyloligellaceae",
  "Labrenzia"               = "Labrenzia",
  "Spiroplasma"             = "Spiroplasma",
  "Blastopirellula"         = "Blastopirellula",
  "SAR116_clade"            = "f_SAR116",
  "Candidatus_Hepatoplasma" = "Candidatus Hepatoplasma",
  "Cyclobacteriaceae"       = "f_Cyclobacteriaceae",
  "Spirochaeta_2"           = "Spirochaeta"
)

# Unified color palette — 20 genera + grey for non-dominant
genus_colors_unified <- c(
  "Endozoicomonas"          = "#1B9E77",
  "c_Gammaproteobacteria"   = "#7570B3",
  "f_Rhodobacteraceae"      = "#E7298A",
  "Oleiphilus"              = "#D95F02",
  "c_Alphaproteobacteria"   = "#377EB8",
  "f_Rhizobiaceae"          = "#66A61E",
  "Mycobacterium"           = "#B2DF8A",
  "Mycoplasma"              = "#E41A1C",
  "c_Bacilli"               = "#CAB2D6",
  "Tistlia"                 = "#6A3D9A",
  "Defluviicoccus"          = "#FF7F00",
  "f_Methyloligellaceae"    = "#FDBF6F",
  "Labrenzia"               = "#B15928",
  "Spiroplasma"             = "#984EA3",
  "Blastopirellula"         = "#A6CEE3",
  "f_SAR116"                = "#33A02C",
  "Candidatus Hepatoplasma" = "#A65628",
  "f_Cyclobacteriaceae"     = "#999999",
  "Spirochaeta"             = "#F781BF",
  "Non-dominant ASVs"       = "grey80"
)

# =============================================================================
# HELPER FUNCTION: Build genus-aggregated data for one species
# =============================================================================

build_species_composition <- function(mt_TSS_obj, species_name, site_levels, site_labels) {
  
  spp_data <- clone(mt_TSS_obj)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == species_name)
  spp_data$sample_table$Year <- as.factor(spp_data$sample_table$Year)
  spp_data$sample_table$Site <- factor(spp_data$sample_table$Site,
                                       levels = site_levels, labels = site_labels)
  spp_data$tidy_dataset()
  
  otu <- spp_data$otu_table
  tax <- spp_data$tax_table
  rel_abund <- sweep(otu, 2, colSums(otu), "/") * 100
  
  asv_genus <- sapply(rownames(tax), function(asv_id) {
    row <- tax[asv_id, ]
    g <- as.character(row[["Genus"]])
    f <- as.character(row[["Family"]])
    cc <- as.character(row[["Class"]])
    
    clean <- function(x, prefix) {
      if (is.na(x) || x == "" || x == prefix) return(NA)
      name <- gsub("^.__", "", x)
      if (grepl("^unclassified|^uncultured|^unknown", name, ignore.case = TRUE))
        return(paste0("uncl. ", gsub("unclassified |uncultured ", "", name)))
      return(name)
    }
    
    label <- clean(g, "g__")
    if (is.na(label)) label <- clean(f, "f__")
    if (is.na(label)) label <- clean(cc, "c__")
    if (is.na(label)) label <- "uncl. Bacteria"
    return(label)
  })
  
  genus_abund <- data.frame(Genus = asv_genus, rel_abund, check.names = FALSE) %>%
    group_by(Genus) %>%
    summarise(across(everything(), sum), .groups = "drop")
  
  meta <- spp_data$sample_table %>% tibble::rownames_to_column("Sample")
  
  long_df <- genus_abund %>%
    pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance") %>%
    left_join(meta, by = "Sample") %>%
    mutate(Species = species_name)
  
  return(long_df)
}

# =============================================================================
# BUILD DATA FOR ALL 4 SPECIES
# =============================================================================

site_levels <- c("BC4", "BC6", "DC1", "DC3")
site_labels <- c("Deerfield Shallow", "Deerfield Deep", "Miami Shallow", "Miami Deep")

all_df <- bind_rows(
  build_species_composition(mt_TSS, "B. asbestinum", site_levels, site_labels),
  build_species_composition(mt_TSS, "E. caribaeorum", site_levels, site_labels),
  build_species_composition(mt_TSS, "E. flexuosa", site_levels, site_labels),
  build_species_composition(mt_TSS, "M. muricata", site_levels, site_labels)
)

# =============================================================================
# APPLY DISPLAY NAMES AND AGGREGATE
# =============================================================================

all_df <- all_df %>%
  mutate(
    Genus_plot = ifelse(Genus %in% top_genera_all,
                        top_genera_display[Genus],
                        "Non-dominant ASVs")
  )

plot_df <- all_df %>%
  group_by(Sample, Genus_plot, Species, Site, Year) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Factor levels
genus_order <- c(rev(unname(top_genera_display)), "Non-dominant ASVs")
plot_df$Genus_plot <- factor(plot_df$Genus_plot, levels = genus_order)

plot_df$Species <- factor(plot_df$Species,
                          levels = c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")
)

# Sort samples by dominant genus within each species/year
dominant_genus <- plot_df %>%
  group_by(Sample) %>%
  slice_max(Abundance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Sample, dominant = Genus_plot, dom_abund = Abundance)

plot_df <- plot_df %>%
  left_join(dominant_genus, by = "Sample") %>%
  arrange(Species, Year, dominant, desc(dom_abund))

sample_order <- plot_df %>%
  distinct(Sample, Species, Year, dominant, dom_abund) %>%
  arrange(Species, Year, dominant, desc(dom_abund)) %>%
  pull(Sample)

plot_df$Sample <- factor(plot_df$Sample, levels = sample_order)

# =============================================================================
# INDIVIDUAL SPECIES PLOTS (no legend, no title) — faceted by Year
# =============================================================================

species_list <- c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")
spp_plots <- list()

for (spp in species_list) {
  
  spp_df <- plot_df %>% filter(Species == spp)
  
  spp_sample_order <- spp_df %>%
    distinct(Sample, Year, dominant, dom_abund) %>%
    arrange(Year, dominant, desc(dom_abund)) %>%
    pull(Sample)
  
  spp_df$Sample <- factor(spp_df$Sample, levels = spp_sample_order)
  
  p <- ggplot(spp_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(~ Year + Site, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = genus_colors_unified, name = NULL, drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.y     = element_blank(),
      axis.text.y      = element_text(size = 10, color = "black"),
      strip.text = element_text(size = 7, face = "bold"),
      strip.background  = element_rect(fill = "white", color = "black"),
      legend.position   = "none",
      panel.spacing     = unit(0.2, "lines"),
      plot.margin       = margin(2, 5, 2, 5)
    )
  
  spp_plots[[spp]] <- p
}

# =============================================================================
# EXTRACT SHARED LEGEND
# =============================================================================

p_legend <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_colors_unified, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
  theme(
    legend.text     = element_text(size = 10, face = "italic"),
    legend.key.size = unit(0.5, "cm")
  )

legend <- get_legend(p_legend)

# =============================================================================
# COMBINE: species labels + plots stacked vertically + shared legend
# =============================================================================

labeled_plots <- list()
for (spp in species_list) {
  label <- ggdraw() + draw_label(spp, fontface = "bold.italic", size = 13, hjust = 0.5)
  labeled_plots[[spp]] <- plot_grid(label, spp_plots[[spp]], ncol = 1, rel_heights = c(0.08, 1))
}

species_column <- plot_grid(
  plotlist = labeled_plots,
  ncol = 1,
  align = "v"
)

p_combined <- plot_grid(
  species_column, legend,
  nrow = 1,
  rel_widths = c(1, 0.25)
)

print(p_combined)

ggsave(
  file.path(output_dir, "Community_composition_ALL_species_byYear.png"),
  p_combined, dpi = 300, width = 16, height = 16, bg = "white"
)
cat("Saved: Community_composition_ALL_species_byYear.png\n")