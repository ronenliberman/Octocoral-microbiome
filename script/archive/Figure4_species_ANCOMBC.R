################################################################################
# Figure 4: Species-Specific Differential Abundance (ANCOM-BC2) — Year Effect
# Author: Ronen Liberman
# Date: March 2026
# Description: 
#   Runs ANCOM-BC2 within each octocoral species comparing years (2024 vs 2023, 
#   2025 vs 2023). Produces a faceted bar plot showing log fold change of 
#   significant ASVs for each species, styled after the ANCOM-BC2 template 
#   (positive/negative LFC colored bars with error bars).
#
# Dependencies: This script is self-contained. It imports QIIME2 artifacts
#   directly and does not depend on other analysis scripts.
################################################################################

# =============================================================================
# SETUP: Packages
# =============================================================================

library(tidyverse)
library(microeco)
library(file2meco)
library(qiime2R)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)
library(ANCOMBC)
library(phyloseq)

set.seed(123)

# =============================================================================
# SETUP: Directories
# =============================================================================

# Working directory — adjust to your machine
setwd("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/")

# Output directory
output_dir <- file.path("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("FIGURE 4: Species-Specific ANCOM-BC2 — Temporal Comparisons\n")
cat("Reference year: 2023\n")
cat("=============================================================================\n\n")

# =============================================================================
# STEP 1: Import QIIME2 data and build microtable
# =============================================================================

cat("--- STEP 1: Importing QIIME2 data ---\n")

data <- qiime2meco(
  "raw data/table-NoEMC.qza",
  taxonomy_table = "raw data/taxonomy.qza",
  sample_table   = "raw data/metadata_all_edited.csv",
  auto_tidy      = TRUE
)

# Filter to prokaryotes
data$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
data$cal_abund()

cat("  Raw data: ", nrow(data$sample_table), "samples,", 
    nrow(data$otu_table), "ASVs\n")

# =============================================================================
# STEP 2: Filter to focal species and sites
# =============================================================================

cat("\n--- STEP 2: Filtering samples ---\n")

d <- clone(data)

# Keep only 4 focal octocoral species (exclude seawater)
d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE"))

# Harmonize site names
d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

# Remove known unreliable samples
samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024")
d$sample_table <- d$sample_table %>%
  filter(!rownames(.) %in% samples_to_remove)

# Rename species codes to full names
species_names <- c(
  "BAST" = "B. asbestinum",
  "ECAR" = "E. caribaeorum",
  "MMUR" = "M. muricata",
  "EFLE" = "E. flexuosa"
)
d$sample_table <- d$sample_table %>%
  mutate(SPP = recode(SPP, !!!species_names))

# Set Year as factor with 2023 as reference level
d$sample_table$Year <- factor(d$sample_table$Year, levels = c("2023", "2024", "2025"))

d$tidy_dataset()

cat("  Filtered data:", nrow(d$sample_table), "samples,",
    nrow(d$otu_table), "ASVs\n")
cat("  Species:", paste(unique(d$sample_table$SPP), collapse = ", "), "\n")
cat("  Years:", paste(levels(d$sample_table$Year), collapse = ", "), "\n\n")

# Check sample counts per species x year
cat("  Sample counts by Species x Year:\n")
print(table(d$sample_table$SPP, d$sample_table$Year))
cat("\n")

# =============================================================================
# STEP 3: Helper — clean taxonomy labels
# =============================================================================

clean_taxonomy <- function(taxa_names) {
  sapply(taxa_names, function(name) {
    parts <- unlist(strsplit(name, "\\|"))
    
    # Try genus first — no prefix
    genus <- parts[grepl("^g__", parts)]
    if (length(genus) > 0) {
      clean_name <- gsub("^g__", "", genus[1])
      if (clean_name != "" && !grepl("unclassified", clean_name, ignore.case = TRUE)) {
        return(clean_name)
      }
    }
    
    # Family — f_ prefix
    family <- parts[grepl("^f__", parts)]
    if (length(family) > 0) {
      clean_name <- gsub("^f__", "", family[1])
      if (clean_name != "" && !grepl("unclassified", clean_name, ignore.case = TRUE)) {
        return(paste0("f_", clean_name))
      }
    }
    
    # Order — o_ prefix
    order_tax <- parts[grepl("^o__", parts)]
    if (length(order_tax) > 0) {
      clean_name <- gsub("^o__", "", order_tax[1])
      if (clean_name != "" && !grepl("unclassified", clean_name, ignore.case = TRUE)) {
        return(paste0("o_", clean_name))
      }
    }
    
    # Class — c_ prefix
    class_tax <- parts[grepl("^c__", parts)]
    if (length(class_tax) > 0) {
      clean_name <- gsub("^c__", "", class_tax[1])
      if (clean_name != "" && !grepl("unclassified", clean_name, ignore.case = TRUE)) {
        return(paste0("c_", clean_name))
      }
    }
    
    # Phylum — p_ prefix
    phylum <- parts[grepl("^p__", parts)]
    if (length(phylum) > 0) {
      clean_name <- gsub("^p__", "", phylum[1])
      if (clean_name != "" && !grepl("unclassified", clean_name, ignore.case = TRUE)) {
        return(paste0("p_", clean_name))
      }
    }
    
    return("Unclassified Bacteria")
  }, USE.NAMES = FALSE)
}

# =============================================================================
# STEP 4: Rename unclassified taxa in taxonomy table
# =============================================================================

rename_unclassified <- function(tax_table_df) {
  tax_df <- as.data.frame(tax_table_df)
  for (i in 1:nrow(tax_df)) {
    unclass_idx <- grep("__$", tax_df[i, ])
    if (length(unclass_idx) > 0) {
      tax_df[i, unclass_idx] <- paste0(
        tax_df[i, unclass_idx],
        "unclassified ", gsub(".__", "", tax_df[i, unclass_idx[1] - 1])
      )
    }
  }
  return(tax_df)
}

# =============================================================================
# STEP 5A: Run ANCOM-BC2 per species — TWO reference levels
#   Run A: Reference = 2023 --> gives 2024 vs 2023, 2025 vs 2023
#   Run B: Reference = 2024 --> gives 2025 vs 2024
# =============================================================================

cat("--- STEP 5: Running ANCOM-BC2 per species (Year effect) ---\n")

species_list <- c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa")

# Storage for both runs
results_ref2023 <- list()
results_ref2024 <- list()

for (spp in species_list) {
  cat("\n  Processing:", spp, "\n")
  
  # --- Subset to species ---
  spp_data <- clone(d)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == spp)
  spp_data$tidy_dataset()
  
  n_samples <- nrow(spp_data$sample_table)
  n_years   <- length(unique(spp_data$sample_table$Year))
  cat("    Samples:", n_samples, "| Years:", n_years, "\n")
  
  if (n_years < 2 || n_samples < 10) {
    cat("    SKIPPED: insufficient samples or years\n")
    next
  }
  
  # Rename unclassified taxa
  spp_data$tax_table <- rename_unclassified(spp_data$tax_table)
  
  # ---- Run A: Reference = 2023 ----
  cat("    ANCOM-BC2 ref=2023...\n")
  spp_data_2023 <- clone(spp_data)
  spp_data_2023$sample_table$Year <- factor(spp_data_2023$sample_table$Year,
                                            levels = c("2023", "2024", "2025"))
  
  tryCatch({
    ancom_2023 <- trans_diff$new(
      dataset         = spp_data_2023,
      method          = "ancombc2",
      group           = "Year",
      taxa_level      = "all",
      filter_thres    = 0.001,
      p_adjust_method = "fdr",
      alpha           = 0.05
    )
    res_2023 <- ancom_2023$res_diff
    res_2023$Species <- spp
    res_2023$ref_year <- "2023"
    n_sig <- sum(res_2023$diff == TRUE & res_2023$Factors != "(Intercept)", na.rm = TRUE)
    cat("      Significant (excl. intercept):", n_sig, "\n")
    results_ref2023[[spp]] <- res_2023
  }, error = function(e) {
    cat("      ERROR:", conditionMessage(e), "\n")
  })
  
  # ---- Run B: Reference = 2024 ----
  cat("    ANCOM-BC2 ref=2024...\n")
  spp_data_2024 <- clone(spp_data)
  spp_data_2024$sample_table$Year <- factor(spp_data_2024$sample_table$Year,
                                            levels = c("2024", "2023", "2025"))
  
  tryCatch({
    ancom_2024 <- trans_diff$new(
      dataset         = spp_data_2024,
      method          = "ancombc2",
      group           = "Year",
      taxa_level      = "all",
      filter_thres    = 0.001,
      p_adjust_method = "fdr",
      alpha           = 0.05
    )
    res_2024 <- ancom_2024$res_diff
    res_2024$Species <- spp
    res_2024$ref_year <- "2024"
    # Keep only 2025 vs 2024 comparison
    res_2024 <- res_2024 %>% filter(Factors %in% c("Year2025", "(Intercept)"))
    n_sig <- sum(res_2024$diff == TRUE & res_2024$Factors != "(Intercept)", na.rm = TRUE)
    cat("      Significant 2025 vs 2024:", n_sig, "\n")
    results_ref2024[[spp]] <- res_2024
  }, error = function(e) {
    cat("      ERROR:", conditionMessage(e), "\n")
  })
}

# =============================================================================
# STEP 5B: Run ANCOM-BC2 per species — Site effect (reference = DC1 / Miami Shallow)
# =============================================================================

cat("\n--- STEP 5B: Running ANCOM-BC2 per species (Site effect) ---\n")
cat("    Reference site: DC1 (Miami Shallow)\n\n")

results_site <- list()

for (spp in species_list) {
  cat("  Processing:", spp, "\n")
  
  spp_data <- clone(d)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == spp)
  
  # Set Site as factor with DC1 (Miami Shallow) as reference
  spp_data$sample_table$Site <- factor(spp_data$sample_table$Site,
                                       levels = c("DC1", "DC3", "BC4", "BC6"))
  spp_data$tidy_dataset()
  
  n_samples <- nrow(spp_data$sample_table)
  n_sites   <- length(unique(spp_data$sample_table$Site))
  cat("    Samples:", n_samples, "| Sites:", n_sites, "\n")
  cat("    Samples per site:\n")
  print(table(spp_data$sample_table$Site))
  
  if (n_sites < 2 || n_samples < 10) {
    cat("    SKIPPED: insufficient samples or sites\n")
    next
  }
  
  # Rename unclassified taxa
  spp_data$tax_table <- rename_unclassified(spp_data$tax_table)
  
  tryCatch({
    ancom_site <- trans_diff$new(
      dataset         = spp_data,
      method          = "ancombc2",
      group           = "Site",
      taxa_level      = "all",
      filter_thres    = 0.001,
      p_adjust_method = "fdr",
      alpha           = 0.05
    )
    
    res_site <- ancom_site$res_diff
    res_site$Species <- spp
    
    n_sig <- sum(res_site$diff == TRUE & res_site$Factors != "(Intercept)", na.rm = TRUE)
    cat("    Significant ASVs (excl. intercept):", n_sig, "\n\n")
    
    results_site[[spp]] <- res_site
    
  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n\n")
  })
}

# =============================================================================
# STEP 6A: Combine and prepare plotting data
# =============================================================================

cat("\n--- STEP 6: Preparing plot data ---\n")

# Combine ref=2023 results (2024 vs 2023, 2025 vs 2023)
combined_2023 <- bind_rows(results_ref2023) %>%
  filter(diff == TRUE, Factors != "(Intercept)") %>%
  mutate(
    Comparison = case_when(
      Factors == "Year2024" ~ "2023 vs 2024",
      Factors == "Year2025" ~ "2023 vs 2025",
      TRUE ~ Factors
    )
  )

# Combine ref=2024 results (2025 vs 2024 only)
combined_2024 <- bind_rows(results_ref2024) %>%
  filter(diff == TRUE, Factors != "(Intercept)") %>%
  mutate(
    Comparison = "2024 vs 2025"
  )

# Merge all significant results
sig_results <- bind_rows(combined_2023, combined_2024) %>%
  mutate(
    clean_tax     = clean_taxonomy(Taxa),
    LFC_direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"),
    ASV_label     = clean_tax
  )

cat("  Total significant ASVs across all comparisons:", nrow(sig_results), "\n")

# Summary per species x comparison
sig_summary <- sig_results %>%
  group_by(Species, Comparison) %>%
  summarise(n_sig = n(), .groups = "drop")
cat("\n  Significant ASVs per species x comparison:\n")
print(sig_summary)

# Save all raw results (both reference levels)
all_raw <- bind_rows(
  bind_rows(results_ref2023),
  bind_rows(results_ref2024)
)
write.csv(all_raw,
          file.path(output_dir, "Figure4_ANCOMBC2_Year_all_results.csv"),
          row.names = FALSE)

# Save significant results
write.csv(sig_results,
          file.path(output_dir, "Figure4_ANCOMBC2_Year_significant.csv"),
          row.names = FALSE)
# =============================================================================
# STEP 6B: Prepare site plotting data
# =============================================================================

cat("\n--- STEP 6B: Preparing site plot data ---\n")

# Site name lookup for clean comparison labels
site_labels <- c(
  "SiteDC3" = "MiD vs MiS",
  "SiteBC4" = "DfS vs MiS",
  "SiteBC6" = "DfD vs MiS"
)

sig_results_site <- bind_rows(results_site) %>%
  filter(diff == TRUE, Factors != "(Intercept)") %>%
  mutate(
    Comparison    = recode(Factors, !!!site_labels, .default = Factors),
    clean_tax     = clean_taxonomy(Taxa),
    LFC_direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"),
    ASV_label     = clean_tax
  )

cat("  Total significant site ASVs:", nrow(sig_results_site), "\n")

sig_summary_site <- sig_results_site %>%
  group_by(Species, Comparison) %>%
  summarise(n_sig = n(), .groups = "drop")
cat("\n  Significant ASVs per species x site comparison:\n")
print(sig_summary_site)

# Save
all_raw_site <- bind_rows(results_site)
write.csv(all_raw_site,
          file.path(output_dir, "Figure4_ANCOMBC2_Site_all_results.csv"),
          row.names = FALSE)
write.csv(sig_results_site,
          file.path(output_dir, "Figure4_ANCOMBC2_Site_significant.csv"),
          row.names = FALSE)

# =============================================================================
# STEP 7 (REVISED): Build GLOBAL ASV numbering across Year AND Site results
# =============================================================================

cat("\n--- STEP 7: Building global ASV numbering and selecting top ASVs ---\n")

n_top <- 15  # max ASVs per species per analysis

# --- 7A: Select top Year ASVs ---
plot_data_year <- sig_results %>%
  group_by(Species) %>%
  arrange(P.adj) %>%
  slice_head(n = n_top) %>%
  ungroup()

# --- 7B: Select top Site ASVs ---
plot_data_site <- sig_results_site %>%
  group_by(Species) %>%
  arrange(P.adj) %>%
  slice_head(n = n_top) %>%
  ungroup()

# --- 7C: Build global numbering for duplicate genera across BOTH panels ---
# Combine all taxa that will appear in either plot
all_plot_taxa <- bind_rows(
  plot_data_year %>% dplyr::select(Species, Taxa, clean_tax),
  plot_data_site %>% dplyr::select(Species, Taxa, clean_tax)
) %>%
  distinct(Species, Taxa, clean_tax)

# For each species, assign consistent numbers to duplicate genera
# Same ASV (same Taxa hash) always gets the same number
global_numbering <- all_plot_taxa %>%
  group_by(Species, clean_tax) %>%
  mutate(
    dup_n    = n(),
    dup_rank = row_number(),
    ASV_label_global = ifelse(dup_n > 1,
                              paste0(clean_tax, " (", dup_rank, ")"),
                              clean_tax)
  ) %>%
  ungroup() %>%
  dplyr::select(Species, Taxa, ASV_label_global)

# Join back to both plot datasets
plot_data_year <- plot_data_year %>%
  left_join(global_numbering, by = c("Species", "Taxa")) %>%
  mutate(ASV_label = ASV_label_global) %>%
  dplyr::select(-ASV_label_global)

plot_data_site <- plot_data_site %>%
  left_join(global_numbering, by = c("Species", "Taxa")) %>%
  mutate(ASV_label = ASV_label_global) %>%
  dplyr::select(-ASV_label_global)

# Handle any remaining duplicates within a single facet cell
# (same label in same Species x Comparison)
deduplicate_within_facet <- function(df) {
  df %>%
    group_by(Species, Comparison, ASV_label) %>%
    mutate(
      facet_dup_n    = n(),
      facet_dup_rank = row_number(),
      ASV_label = ifelse(facet_dup_n > 1,
                         paste0(ASV_label, ".", facet_dup_rank),
                         ASV_label)
    ) %>%
    ungroup() %>%
    dplyr::select(-facet_dup_n, -facet_dup_rank)
}

plot_data_year <- deduplicate_within_facet(plot_data_year)
plot_data_site <- deduplicate_within_facet(plot_data_site)

cat("  Year ASVs for plot:", nrow(plot_data_year), "\n")
cat("  Site ASVs for plot:", nrow(plot_data_site), "\n\n")

# =============================================================================
# STEP 8A: Generate Year panel (Panel A)
# =============================================================================

cat("--- STEP 8A: Generating Year panel (Panel A) ---\n")

lfc_colors <- c("Positive LFC" = "#FA8072",
                "Negative LFC" = "#66CDAA")

# Factor levels
plot_data_year$Species <- factor(plot_data_year$Species,
                                 levels = c("B. asbestinum", "E. caribaeorum",
                                            "E. flexuosa", "M. muricata"))

# Chronological order: early year first
plot_data_year$Comparison <- factor(plot_data_year$Comparison,
                                    levels = c("2023 vs 2024", "2023 vs 2025", "2024 vs 2025"))

# Facet-specific ordering
plot_data_year <- plot_data_year %>%
  mutate(facet_label = paste0(Species, "___", ASV_label),
         facet_label = reorder(facet_label, lfc))

panel_A <- ggplot(plot_data_year, aes(x = facet_label, y = lfc, fill = LFC_direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se),
                width = 0.25, linewidth = 0.4, color = "black") +
  facet_grid(Species ~ Comparison, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) gsub("^.*___", "", x)) +
  scale_fill_manual(values = lfc_colors, name = "") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  labs(x = NULL, y = "Log fold change") +
  theme_classic(base_size = 12) +
  theme(
    strip.text.y = element_text(size = 15, face = "bold.italic", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, face = "italic"),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.6, "lines"),
    plot.clip = "off",
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    plot.margin = margin(5, 10, 20, 5)
  )

cat("  Panel A created.\n")

# =============================================================================
# STEP 8B: Generate Site panel (Panel B)
# =============================================================================

cat("--- STEP 8B: Generating Site panel (Panel B) ---\n")

plot_data_site$Species <- factor(plot_data_site$Species,
                                 levels = c("B. asbestinum", "E. caribaeorum",
                                            "E. flexuosa", "M. muricata"))
plot_data_site$Comparison <- factor(plot_data_site$Comparison,
                                    levels = c("MiD vs MiS",
                                               "DfS vs MiS",
                                               "DfD vs MiS"))
plot_data_site$Species <- droplevels(plot_data_site$Species)
plot_data_site <- plot_data_site %>%
  mutate(facet_label = paste0(Species, "___", ASV_label),
         facet_label = reorder(facet_label, lfc))

panel_B <- ggplot(plot_data_site, aes(x = facet_label, y = lfc, fill = LFC_direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se),
                width = 0.25, linewidth = 0.4, color = "black") +
  facet_grid(Species ~ Comparison, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) gsub("^.*___", "", x)) +
  scale_fill_manual(values = lfc_colors, name = "") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  theme_classic(base_size = 12) +
  theme(
    strip.text.y = element_text(size = 15, face = "bold.italic", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, face = "italic"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.6, "lines"),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    plot.margin = margin(5, 10, 5, 5)
  )

cat("  Panel B created.\n")

# =============================================================================
# STEP 9: Combine into one figure with shared legend
# =============================================================================

cat("\n--- STEP 9: Combining panels ---\n")

# Extract legend from a temporary plot
legend_plot <- ggplot(plot_data_year, aes(x = facet_label, y = lfc, fill = LFC_direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = lfc_colors, name = "") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.72, "cm"),
        legend.spacing.x = unit(0.5, "cm"))

shared_legend <- cowplot::get_legend(legend_plot)

# Horizontal separator line
separator <- ggplot() +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0,
           color = "grey50", linewidth = 0.5) +
  theme_void() +
  theme(plot.margin = margin(2, 0, 2, 0))

# Stack: Panel A, separator, Panel B, legend at bottom
fig4_combined <- plot_grid(
  panel_A,
  separator,
  panel_B,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.02, 1, 0.06),
  align = "v",
  axis = "lr"
)

fig4_combined <- ggdraw(fig4_combined) +
  draw_label("A)", x = 0.01, y = 0.99, fontface = "bold", size = 16, hjust = 0) +
  draw_label("B)", x = 0.01, y = 0.5, fontface = "bold", size = 16, hjust = 0)
print(fig4_combined)

# Dynamic sizing
n_year_bars <- nrow(plot_data_year)
n_site_bars <- nrow(plot_data_site)
total_bars  <- n_year_bars + n_site_bars
plot_height <- max(14, total_bars * 0.3 + 6)
plot_width  <- max(18, max(
  length(unique(plot_data_year$Comparison)),
  length(unique(plot_data_site$Comparison))
) * 6 + 4)

ggsave(file.path(output_dir, "Figure4_ANCOMBC2_combined.png"),
       fig4_combined, dpi = 300, width = 16, height = 22, bg = "white")
ggsave(file.path(output_dir, "Figure4_ANCOMBC2_combined.pdf"),
       fig4_combined, width = 16, height = 22, bg = "white")

cat("  Combined figure saved.\n")

# =============================================================================
# STEP 10: Summary statistics and CSV export
# =============================================================================

cat("\n--- STEP 10: Summary ---\n\n")

# Year summary
summary_year <- sig_results %>%
  group_by(Species, Comparison, LFC_direction) %>%
  summarise(n_ASVs = n(), mean_LFC = round(mean(lfc), 2),
            median_LFC = round(median(lfc), 2), .groups = "drop") %>%
  mutate(Analysis = "Temporal") %>%
  arrange(Species, Comparison)

# Site summary
summary_site <- sig_results_site %>%
  group_by(Species, Comparison, LFC_direction) %>%
  summarise(n_ASVs = n(), mean_LFC = round(mean(lfc), 2),
            median_LFC = round(median(lfc), 2), .groups = "drop") %>%
  mutate(Analysis = "Spatial") %>%
  arrange(Species, Comparison)

summary_all <- bind_rows(summary_year, summary_site)
cat("Combined summary of significant differentially abundant ASVs:\n")
print(summary_all)

write.csv(summary_all,
          file.path(output_dir, "Figure4_ANCOMBC2_combined_summary_stats.csv"),
          row.names = FALSE)

cat("\n=============================================================================\n")
cat("FIGURE 4 COMPLETE\n")
cat("=============================================================================\n")
cat("\nOutput files:\n")
cat("  - Figure4_ANCOMBC2_Year_all_results.csv\n")
cat("  - Figure4_ANCOMBC2_Year_significant.csv\n")
cat("  - Figure4_ANCOMBC2_Site_all_results.csv\n")
cat("  - Figure4_ANCOMBC2_Site_significant.csv\n")
cat("  - Figure4_ANCOMBC2_combined_summary_stats.csv\n")
cat("  - Figure4_ANCOMBC2_combined.png/pdf\n")