# =============================================================================
# REVISED CORE MICROBIOME VISUALIZATION - ALL THRESHOLDS
# Addresses reviewer comments:
#   - Show all threshold results (75%, 90%, 100%), not just 90%
#   - Report at ASV level consistently
#   - Heatmap-style figure inspired by seed lot/performance figure
# =============================================================================

# --- PACKAGES ----------------------------------------------------------------
library(microeco)
library(qiime2R)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(writexl)
library(cowplot)
library(grid)
library(gridExtra)

set.seed(123)

# --- Adjust this path to your working directory ---
setwd("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/")
# --- OUTPUT DIRECTORY --------------------------------------------------------
output_dir <- "~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# SECTION 1: DATA IMPORT & SETUP (qiime2R → microeco)
# =============================================================================

# Import QIIME2 data
data <- qiime2meco(
  "raw_data/table-NoEMC-PEE2.qza",
  taxonomy_table = "raw_data/taxonomy-PEE2.qza",
  sample_table   = "raw_data/metadata_all_edited.csv",
  auto_tidy      = TRUE
)

# Filter to prokaryotes
data$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
data$cal_abund()

# Clone and filter to focal species + sites
d <- clone(data)
d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE", "SW"))
d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

# Remove unreliable samples
samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024", "EFLE_DC3_3_2024")
d$sample_table <- d$sample_table %>%
  filter(!rownames(.) %in% samples_to_remove)

# Rename species
species.names <- list(
  "SW"   = "Seawater",
  "BAST" = "B. asbestinum",
  "ECAR" = "E. caribaeorum",
  "MMUR" = "M. muricata",
  "EFLE" = "E. flexuosa"
)
d$sample_table <- d$sample_table %>%
  mutate(SPP = recode(SPP, !!!species.names))
d$tidy_dataset()

# --- TSS normalization (for relative abundance heatmap) ----------------------
tmp_tss <- trans_norm$new(dataset = d)
mt_tss  <- tmp_tss$norm(method = "TSS")

# --- Rarefaction (for prevalence-based core identification) ------------------
tmp_rare   <- trans_norm$new(dataset = d)
mt_rarefied <- tmp_rare$norm(method = "rarefy", sample.size = 20000)

cat("TSS dataset:", ncol(mt_tss$otu_table), "samples,", nrow(mt_tss$otu_table), "ASVs\n")
cat("Rarefied dataset:", ncol(mt_rarefied$otu_table), "samples,", nrow(mt_rarefied$otu_table), "ASVs\n")

# =============================================================================
# SECTION 2: CORE MICROBIOME CALCULATION AT ALL THRESHOLDS
# =============================================================================

target_species        <- c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa")
prevalence_thresholds <- c(80, 90, 100)  # Changed 75 → 80

# Function to calculate prevalence-based core
calculate_prevalence_core <- function(microtable_obj, prev_threshold) {
  otu_table <- microtable_obj$otu_table
  tax_table <- microtable_obj$tax_table
  
  presence_absence <- otu_table
  presence_absence[presence_absence > 0] <- 1
  prevalence <- rowSums(presence_absence) / ncol(presence_absence) * 100
  
  core_taxa <- names(prevalence[prevalence >= prev_threshold])
  
  if (length(core_taxa) > 0) {
    core_otu <- otu_table[core_taxa, , drop = FALSE]
    total_reads_per_sample <- colSums(otu_table)
    core_rel_abund <- sweep(core_otu, 2, total_reads_per_sample, "/") * 100
    mean_rel_abundance <- rowMeans(core_rel_abund)
    
    core_summary <- data.frame(
      ASV_ID = core_taxa,
      Prevalence = prevalence[core_taxa],
      Mean_Relative_Abundance_Percent = mean_rel_abundance,
      tax_table[core_taxa, ],
      stringsAsFactors = FALSE
    )
    core_summary <- core_summary[order(-core_summary$Mean_Relative_Abundance_Percent), ]
  } else {
    core_summary <- data.frame()
  }
  
  return(list(
    core_taxa    = core_taxa,
    core_summary = core_summary,
    n_core_taxa  = length(core_taxa)
  ))
}

# Fix unclassified taxonomy labels in rarefied data
core_data <- clone(mt_tss)
rename_unclassified <- core_data$tax_table
for (i in 1:nrow(rename_unclassified)) {
  if (any(grepl("__$", rename_unclassified[i, ]))) {
    matchall <- grep("__$", rename_unclassified[i, ])
    rename_unclassified[i, matchall] <- paste0(
      rename_unclassified[i, matchall],
      "unclassified ",
      gsub(".__", "", rename_unclassified[i, matchall[1] - 1])
    )
  }
}
core_data$tax_table <- rename_unclassified

# Calculate core at all thresholds for all species
all_core_results <- list()

for (threshold in prevalence_thresholds) {
  cat("\n--- Prevalence threshold:", threshold, "% ---\n")
  threshold_results <- list()
  
  for (species in target_species) {
    species_data <- clone(core_data)
    species_data$sample_table %<>% subset(SPP == species)
    species_data$tidy_dataset()
    
    core_result <- calculate_prevalence_core(species_data, threshold)
    threshold_results[[species]] <- core_result
    cat(species, ":", core_result$n_core_taxa, "core ASVs\n")
  }
  all_core_results[[as.character(threshold)]] <- threshold_results
}

# =============================================================================
# SECTION 3: BUILD COMBINED DATA FRAME FOR ALL CORE ASVs × SPECIES × THRESHOLDS
# =============================================================================

combined_core_df <- data.frame()

for (threshold in prevalence_thresholds) {
  for (species in target_species) {
    core_result <- all_core_results[[as.character(threshold)]][[species]]
    if (core_result$n_core_taxa > 0) {
      tmp_df <- core_result$core_summary %>%
        mutate(
          Species   = species,
          Threshold = threshold
        )
      combined_core_df <- bind_rows(combined_core_df, tmp_df)
    }
  }
}

# Convert Rle columns
combined_core_df[] <- lapply(combined_core_df, function(x) if (is(x, "Rle")) as.vector(x) else x)
combined_core_df <- as.data.frame(combined_core_df)

# Clean genus names for display
combined_core_df <- combined_core_df %>%
  mutate(
    Genus_clean = gsub("^g__", "", Genus),
    Family_clean = gsub("^f__", "", Family),
    Display_Label = case_when(
      grepl("unclassified", Genus_clean) & !grepl("unclassified", Family_clean) ~
        paste0("Uncl. ", Family_clean),
      grepl("unclassified", Genus_clean) & grepl("unclassified", Family_clean) ~
        paste0("Uncl. ", gsub("^c__", "", Class)),
      TRUE ~ Genus_clean
    )
  )

# Create unique ASV labels
combined_core_df <- combined_core_df %>%
  group_by(Species, Display_Label) %>%
  mutate(
    n_asvs_same_genus = n_distinct(ASV_ID),
    asv_rank = dense_rank(ASV_ID)
  ) %>%
  ungroup() %>%
  mutate(
    ASV_Label = ifelse(n_asvs_same_genus > 1,
                       paste0(Display_Label, " ", asv_rank),
                       Display_Label)
  )

# Threshold membership
threshold_membership <- as.data.frame(combined_core_df) %>%
  select(ASV_ID, Species, Threshold) %>%
  distinct() %>%
  group_by(ASV_ID, Species) %>%
  summarise(
    Min_Threshold = min(Threshold),
    Max_Threshold = max(Threshold),
    Strictest_Threshold = max(Threshold),
    Threshold_Label = paste0(sort(unique(Threshold)), "%", collapse = " + "),
    .groups = "drop"
  )

# Plot data
plot_data <- as.data.frame(combined_core_df) %>%
  group_by(ASV_ID, Species) %>%
  filter(Threshold == max(Threshold)) %>%
  ungroup() %>%
  distinct(ASV_ID, Species, .keep_all = TRUE) %>%
  left_join(threshold_membership, by = c("ASV_ID", "Species"))

cat("\nTotal unique ASV x Species combinations:", nrow(plot_data), "\n")
cat("Unique ASVs across all species:", n_distinct(plot_data$ASV_ID), "\n")

# =============================================================================
# SECTION 4: PREPARE DATA FOR HEATMAP
# =============================================================================

all_core_asvs <- unique(plot_data$ASV_ID)

tss_coral <- clone(mt_tss)
tss_coral$sample_table %<>% subset(SPP %in% target_species)
tss_coral$tidy_dataset()

core_rel_abund_matrix <- tss_coral$otu_table[
  rownames(tss_coral$otu_table) %in% all_core_asvs, , drop = FALSE
] * 100

cat("Core ASVs found in TSS data:", nrow(core_rel_abund_matrix), "of", length(all_core_asvs), "\n")

heatmap_long <- core_rel_abund_matrix %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  pivot_longer(-ASV_ID, names_to = "SampleID", values_to = "Rel_Abundance") %>%
  left_join(
    tss_coral$sample_table %>%
      as.data.frame() %>%
      rownames_to_column("SampleID") %>%
      select(SampleID, SPP, Year, Site),
    by = "SampleID"
  )

# =============================================================================
# SECTION 5: HEATMAP-STYLE VISUALIZATION
# =============================================================================

species_colors <- c(
  "B. asbestinum"  = "#E31A1C",
  "E. caribaeorum" = "#FF7F00",
  "M. muricata"    = "#984EA3",
  "E. flexuosa"    = "#4DAF4A"
)

species_order_levels <- c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa")

# --- Identify shared ASVs (core in 2+ species at ANY threshold) ---
all_core_asv_species <- data.frame()
for (thresh in prevalence_thresholds) {
  for (sp in target_species) {
    core_asvs <- all_core_results[[as.character(thresh)]][[sp]]$core_taxa
    if (length(core_asvs) > 0) {
      all_core_asv_species <- bind_rows(all_core_asv_species, data.frame(
        ASV_ID = core_asvs,
        Species = sp,
        Threshold = thresh,
        stringsAsFactors = FALSE
      ))
    }
  }
}

shared_asvs <- all_core_asv_species %>%
  select(ASV_ID, Species) %>%
  distinct() %>%
  group_by(ASV_ID) %>%
  summarise(
    n_species = n_distinct(Species),
    shared_with = paste(sort(unique(Species)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_species >= 2)

cat("\n=== SHARED CORE ASVs (present in 2+ species) ===\n")
cat("Number of shared ASVs:", nrow(shared_asvs), "\n")
print(shared_asvs)

# --- Prepare per-species ASV ordering ---
species_asv_order <- plot_data %>%
  arrange(Species, desc(Mean_Relative_Abundance_Percent)) %>%
  group_by(Species) %>%
  mutate(ASV_Order = row_number()) %>%
  ungroup()

# Y_Label: just the ASV label, no hash, no species tag
species_asv_order <- species_asv_order %>%
  mutate(Y_Label = ASV_Label)

# Add shared ASV flag with * / ** distinction
species_asv_order <- species_asv_order %>%
  left_join(shared_asvs %>% select(ASV_ID, n_species, shared_with), by = "ASV_ID") %>%
  mutate(
    is_shared = !is.na(n_species) & n_species >= 2,
    n_species = ifelse(is.na(n_species), 1, n_species),
    shared_symbol = case_when(
      n_species >= 3 ~ "**",
      n_species == 2 ~ "*",
      TRUE ~ ""
    )
  )

# Global y-axis order (used only as reference; droplevels handles per-species)
y_label_order_vec <- species_asv_order %>%
  mutate(Species = factor(Species, levels = species_order_levels)) %>%
  arrange(Species, ASV_Order) %>%
  pull(Y_Label)

# Join heatmap data with ordering
heatmap_plot_data <- heatmap_long %>%
  inner_join(
    species_asv_order %>%
      select(ASV_ID, Species, ASV_Order, Y_Label, ASV_Label,
             Strictest_Threshold, Threshold_Label,
             Mean_Relative_Abundance_Percent, Prevalence,
             is_shared, n_species, shared_symbol),
    by = c("ASV_ID", "SPP" = "Species")
  )

heatmap_plot_data <- heatmap_plot_data %>%
  arrange(SPP, Year, Site) %>%
  group_by(SPP) %>%
  mutate(Sample_Order = as.numeric(factor(SampleID, levels = unique(SampleID)))) %>%
  ungroup()

# =============================================================================
# SECTION 6: BUILD THE COMPOSITE FIGURE
# =============================================================================

# ASV counts per species per threshold for legend
asv_counts_per_species <- data.frame()
for (sp in species_order_levels) {
  counts <- sapply(prevalence_thresholds, function(thresh) {
    all_core_results[[as.character(thresh)]][[sp]]$n_core_taxa
  })
  asv_counts_per_species <- bind_rows(asv_counts_per_species, data.frame(
    Species = sp,
    label = paste0(sp, " (", paste(counts, collapse = ", "), ")"),
    stringsAsFactors = FALSE
  ))
}

build_species_panel <- function(sp_name, sp_data, sp_asv_info, is_last = FALSE) {
  
  # Order Y_Label by ASV_Order within this species, convert to factor
  label_order <- sp_asv_info %>%
    arrange(ASV_Order) %>%
    pull(Y_Label) %>%
    unique()
  
  sp_data$Y_Label <- factor(sp_data$Y_Label, levels = rev(label_order))
  sp_asv_info$Y_Label <- factor(sp_asv_info$Y_Label, levels = rev(label_order))
  
  # --- Panel A: Heatmap ---
  p_heat <- ggplot(sp_data, aes(x = reorder(SampleID, Sample_Order),
                                y = Y_Label,
                                fill = Rel_Abundance)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(
      colors = c("white", "#FFF5EB", "#FDD49E", "#FDBB84",
                 "#FC8D59", "#EF6548", "#D7301F", "#990000"),
      values = scales::rescale(c(0, 0.5, 1, 2, 5, 10, 20, 40)),
      limits = c(0, NA),
      na.value = "white",
      name = "Relative\nAbundance (%)"
    ) +
    labs(x = NULL, y = NULL, title = sp_name) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10, face = "italic", color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 13, face = "bold.italic", hjust = 0),
      legend.position = "none",
      plot.margin = margin(2, 0, 2, 2)
    )
  
  # --- Panel B: Threshold indicator tiles ---
  asv_ids <- sp_asv_info$ASV_ID
  
  real_threshold_data <- data.frame()
  for (thresh in c(80, 90, 100)) {
    core_asvs_at_thresh <- all_core_results[[as.character(thresh)]][[sp_name]]$core_taxa
    for (asv in asv_ids) {
      real_threshold_data <- bind_rows(real_threshold_data, data.frame(
        ASV_ID = asv,
        Threshold_level = paste0(thresh, "%"),
        Present = asv %in% core_asvs_at_thresh,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  real_threshold_data <- real_threshold_data %>%
    left_join(sp_asv_info %>% select(ASV_ID, Y_Label), by = "ASV_ID") %>%
    mutate(Threshold_level = factor(Threshold_level, levels = c("80%", "90%", "100%")))
  
  p_thresh <- ggplot(real_threshold_data, aes(x = Threshold_level, y = Y_Label)) +
    geom_point(aes(color = Threshold_level, size = Present), shape = 15) +
    scale_color_manual(
      values = c("80%" = "#A6CEE3", "90%" = "#1F78B4", "100%" = "#08306B"),
      name = "Threshold"
    ) +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 0), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 0)
    )
  
  # --- Panel C: Mean relative abundance bars + shared ASV symbols ---
  max_abund <- max(sp_asv_info$Mean_Relative_Abundance_Percent, na.rm = TRUE)
  x_upper <- max(max_abund * 1.2, 1)
  
  p_abund <- ggplot(sp_asv_info, aes(x = Mean_Relative_Abundance_Percent, y = Y_Label)) +
    geom_col(fill = species_colors[sp_name], alpha = 0.7, width = 0.7) +
    geom_text(
      data = sp_asv_info %>% filter(is_shared),
      aes(x = Mean_Relative_Abundance_Percent + (max_abund * 0.05),
          y = Y_Label, label = shared_symbol),
      size = 5, fontface = "bold", hjust = 0, vjust = 0.75
    ) +
    scale_x_continuous(limits = c(0, x_upper), expand = c(0, 0)) +
    labs(
      x = if (is_last) "Mean Rel. Abundance (%)" else NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = if (is_last) element_text(size = 11, face = "bold", color = "black") else element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 0)
    )
  
  list(heat = p_heat, thresh = p_thresh, abund = p_abund)
}

# Build panels — last species gets the x-axis title
all_panels <- list()
for (i in seq_along(species_order_levels)) {
  sp <- species_order_levels[i]
  is_last <- (i == length(species_order_levels))
  
  sp_heatmap_data <- heatmap_plot_data %>% filter(SPP == sp)
  sp_asv_info <- species_asv_order %>% filter(Species == sp)
  
  if (nrow(sp_heatmap_data) > 0) {
    all_panels[[sp]] <- build_species_panel(sp, sp_heatmap_data, sp_asv_info, is_last)
  }
}

# =============================================================================
# SECTION 7: ASSEMBLE FINAL COMPOSITE FIGURE
# =============================================================================

build_species_row <- function(panels) {
  plot_grid(
    panels$heat, panels$thresh, panels$abund,
    nrow = 1,
    rel_widths = c(0.7, 0.12, 0.18),
    align = "h", axis = "tb"
  )
}

species_rows <- lapply(all_panels, build_species_row)

n_asvs_per_species <- species_asv_order %>%
  group_by(Species) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(match(Species, species_order_levels))

rel_heights <- pmax(n_asvs_per_species$n + 2, 5)

p_combined <- plot_grid(
  plotlist = species_rows,
  ncol = 1,
  rel_heights = rel_heights,
  align = "v"
)

# --- Legends (bottom-right, vertically aligned) ---

p_legend_heat <- ggplot(data.frame(x = 1:8, fill = c(0, 0.5, 1, 2, 5, 10, 20, 40)),
                        aes(x, x, fill = fill)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "#FFF5EB", "#FDD49E", "#FDBB84",
               "#FC8D59", "#EF6548", "#D7301F", "#990000"),
    values = scales::rescale(c(0, 0.5, 1, 2, 5, 10, 20, 40)),
    name = "Relative\nAbundance (%)",
    breaks = c(0, 5, 10, 20, 40)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 17)
  )

p_legend_thresh <- ggplot(data.frame(x = 1:3,
                                     thresh = factor(c("80%", "90%", "100%"),
                                                     levels = c("80%", "90%", "100%"))),
                          aes(x, 1, color = thresh)) +
  geom_point(shape = 15, size = 5) +
  scale_color_manual(
    values = c("80%" = "#A6CEE3", "90%" = "#1F78B4", "100%" = "#08306B"),
    name = "Prevalence\nThreshold"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 17)
  )

species_legend_labels_expr <- setNames(
  asv_counts_per_species$label,
  asv_counts_per_species$Species
)

p_legend_sp <- ggplot(data.frame(x = 1:4,
                                 sp = factor(species_order_levels, levels = species_order_levels)),
                      aes(x, 1, fill = sp)) +
  geom_tile() +
  scale_fill_manual(
    values = species_colors,
    name = "Species",
    labels = species_legend_labels_expr
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "italic")
  )

leg_heat   <- get_legend(p_legend_heat)
leg_thresh <- get_legend(p_legend_thresh)
leg_sp     <- get_legend(p_legend_sp)

# Stack legends vertically, pushed to bottom with NULL spacer on top
legends_combined <- plot_grid(
  NULL,
  leg_heat, leg_thresh, leg_sp,
  ncol = 1,
  rel_heights = c(1, 0.4, 0.25, 0.35),
  align = "v"
)

# Final figure
p_final <- plot_grid(
  p_combined, legends_combined,
  nrow = 1, rel_widths = c(0.82, 0.18)
)

# Save
ggsave(file.path(output_dir, "Core_Microbiome_AllThresholds_Heatmap_TSS.png"),
       p_final, width = 16, height = 20, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "Core_Microbiome_AllThresholds_Heatmap_TSS.pdf"),
       p_final, width = 16, height = 14)

cat("\n=== Updated figure saved ===\n")

# =============================================================================
# SECTION 8: SUMMARY TABLE FOR SUPPLEMENT
# =============================================================================

summary_table <- data.frame()
for (sp in target_species) {
  for (thresh in prevalence_thresholds) {
    core_result <- all_core_results[[as.character(thresh)]][[sp]]
    if (core_result$n_core_taxa > 0) {
      tmp <- core_result$core_summary %>%
        as.data.frame() %>%
        mutate(Species = sp, Threshold = thresh)
      summary_table <- bind_rows(summary_table, tmp)
    }
  }
}

threshold_wide <- summary_table %>%
  select(ASV_ID, Species, Threshold) %>%
  distinct() %>%
  mutate(Core = TRUE) %>%
  pivot_wider(names_from = Threshold, values_from = Core,
              names_prefix = "Core_", values_fill = FALSE)

best_info <- summary_table %>%
  as.data.frame() %>%
  group_by(ASV_ID, Species) %>%
  filter(Threshold == max(Threshold)) %>%
  ungroup() %>%
  distinct(ASV_ID, Species, .keep_all = TRUE)

final_table <- best_info %>%
  left_join(threshold_wide, by = c("ASV_ID", "Species")) %>%
  left_join(shared_asvs %>% select(ASV_ID, n_species, shared_with), by = "ASV_ID") %>%
  mutate(
    Shared_across_species = ifelse(!is.na(n_species) & n_species >= 2, "Yes", "No"),
    Shared_with = ifelse(is.na(shared_with), "", shared_with)
  ) %>%
  select(-n_species, -shared_with) %>%
  arrange(Species, desc(Mean_Relative_Abundance_Percent))

write.csv(final_table,
          file.path(output_dir, "Core_Microbiome_AllThresholds_Summary_TSS.csv"),
          row.names = FALSE)

cat("\nSummary table saved: Core_Microbiome_AllThresholds_Summary_TSS.csv\n")

# =============================================================================
# SECTION 9: PRINT SUMMARY STATISTICS
# =============================================================================

cat("\n=== CORE MICROBIOME SUMMARY ACROSS ALL THRESHOLDS ===\n\n")

for (sp in target_species) {
  cat(sp, ":\n")
  for (thresh in prevalence_thresholds) {
    n <- all_core_results[[as.character(thresh)]][[sp]]$n_core_taxa
    cat("  ", thresh, "% threshold:", n, "core ASVs\n")
  }
  cat("\n")
}

cat("\n=== SHARED CORE ASVs ===\n")
cat("ASVs found in the core microbiome of 2+ species:\n\n")
for (i in 1:nrow(shared_asvs)) {
  asv <- shared_asvs$ASV_ID[i]
  tax_info <- as.data.frame(combined_core_df) %>%
    filter(ASV_ID == asv) %>%
    distinct(ASV_ID, .keep_all = TRUE)
  symbol <- ifelse(shared_asvs$n_species[i] >= 3, "**", "*")
  cat("  ", symbol, " ", substr(asv, 1, 8), "... -",
      tax_info$Display_Label[1],
      "- shared by:", shared_asvs$shared_with[i], "\n")
}