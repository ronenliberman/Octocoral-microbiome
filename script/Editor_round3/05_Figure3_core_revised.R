# =============================================================================
# 05_Figure3_core_revised.R — Figure 3 (core microbiome) without empty rows
# =============================================================================
# Editor comment addressed: "Fig 3 has too much empty space. Please only show
# ASVs representing some threshold such as >1% of the community."
#
# Changes from script/Core_Microbiome_Revised_AllThresholds.R (R1):
#   1. Core membership is computed exactly as before (80/90/100% prevalence,
#      TSS data), but only core ASVs with mean relative abundance >=
#      MIN_MEAN_RELABUND_FIG3 (1%) within the host are drawn. The complete
#      core list stays in Appendix 2 (written again here).
#   2. Species colours fixed (R1 swapped M. muricata and E. flexuosa).
#   3. Row labels come from ASV_stable_IDs.csv (02), so an ASV has the same
#      name here as in Figs 4–6 and the appendices.
#   4. Legend counts still give the FULL number of core ASVs at 80/90/100%;
#      the number of rows drawn is reported separately.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

suppressPackageStartupMessages({
  library(cowplot)
})

cat("\n========== 05: Figure 3 (core microbiome, >= 1% rows) ==========\n\n")

# =============================================================================
# SECTION 1: Data and TSS normalisation
# =============================================================================

d <- load_octocoral_dataset(include_seawater = FALSE)
mt_tss <- trans_norm$new(dataset = d)$norm(method = "TSS")

ids <- read.csv(file.path(out_dir, "ASV_stable_IDs.csv"), stringsAsFactors = FALSE)

target_species        <- c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa")
prevalence_thresholds <- c(80, 90, 100)

# =============================================================================
# SECTION 2: Core ASVs per species and threshold (unchanged method)
# =============================================================================

calculate_prevalence_core <- function(mt, prev_threshold) {
  otu <- mt$otu_table
  prevalence <- rowSums(otu > 0) / ncol(otu) * 100
  core <- names(prevalence[prevalence >= prev_threshold])
  rel <- sweep(otu[core, , drop = FALSE], 2, colSums(otu), "/") * 100
  tibble(ASV_ID = core,
         Prevalence = prevalence[core],
         Mean_Relative_Abundance_Percent = rowMeans(rel))
}

core_long <- bind_rows(lapply(target_species, function(sp) {
  mt <- clone(mt_tss)
  mt$sample_table %<>% subset(SPP == sp)
  mt$tidy_dataset()
  bind_rows(lapply(prevalence_thresholds, function(th) {
    calculate_prevalence_core(mt, th) %>% mutate(Species = sp, Threshold = th)
  }))
}))

core_counts <- core_long %>% count(Species, Threshold)
print(core_counts %>% pivot_wider(names_from = Threshold, values_from = n))

# One row per ASV x species, with the thresholds it meets
core_asv <- core_long %>%
  group_by(Species, ASV_ID) %>%
  summarise(Prevalence = dplyr::first(Prevalence),
            Mean_Relative_Abundance_Percent = dplyr::first(Mean_Relative_Abundance_Percent),
            Core_80 = 80 %in% Threshold, Core_90 = 90 %in% Threshold,
            Core_100 = 100 %in% Threshold, .groups = "drop")

# Shared = core (any threshold) in 2+ host species
shared <- core_asv %>% count(ASV_ID, name = "n_species")
core_asv <- core_asv %>%
  left_join(shared, by = "ASV_ID") %>%
  mutate(shared_symbol = case_when(n_species >= 3 ~ "**", n_species == 2 ~ "*", TRUE ~ ""))

# Labels: stable IDs, falling back to lowest SILVA rank for any ASV missing from 02
# Only the rank columns: tax_table also has a "Species" (SILVA species) column,
# which would clash with the host-species column of the same name.
tax_df <- d$tax_table %>% as.data.frame() %>% rownames_to_column("ASV_ID") %>%
  dplyr::select(ASV_ID, Phylum, Class, Order, Family, Genus)
core_asv <- core_asv %>%
  left_join(ids %>% dplyr::select(ASV_ID, Figure_Label), by = "ASV_ID") %>%
  left_join(tax_df, by = "ASV_ID") %>%
  mutate(Figure_Label = coalesce(Figure_Label, paste0(strip_prefix(Genus), " [", substr(ASV_ID, 1, 6), "]")))

# Appendix 2 (complete core list, every threshold) — nothing filtered
appendix2 <- core_asv %>%
  dplyr::select(Species, ASV_ID, Figure_Label, Prevalence, Mean_Relative_Abundance_Percent,
                Core_80, Core_90, Core_100, n_species, Phylum, Class, Order, Family, Genus) %>%
  mutate(Shown_in_Fig3 = Mean_Relative_Abundance_Percent >= MIN_MEAN_RELABUND_FIG3) %>%
  arrange(factor(Species, levels = target_species), desc(Mean_Relative_Abundance_Percent))
write.csv(appendix2, file.path(out_dir, "Appendix2_core_microbiome_all_thresholds.csv"), row.names = FALSE)

# =============================================================================
# SECTION 3: Rows to draw (>= 1% mean relative abundance)
# =============================================================================

plot_asv <- core_asv %>%
  filter(Mean_Relative_Abundance_Percent >= MIN_MEAN_RELABUND_FIG3) %>%
  group_by(Species) %>%
  arrange(desc(Mean_Relative_Abundance_Percent), .by_group = TRUE) %>%
  mutate(ASV_Order = row_number()) %>%
  ungroup()

cat("\nRows drawn per species (>=", MIN_MEAN_RELABUND_FIG3, "%):\n")
print(plot_asv %>% count(Species))

tss_coral <- clone(mt_tss)
tss_coral$sample_table %<>% subset(SPP %in% target_species)
tss_coral$tidy_dataset()

heat_long <- (tss_coral$otu_table[unique(plot_asv$ASV_ID), , drop = FALSE] * 100) %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  pivot_longer(-ASV_ID, names_to = "SampleID", values_to = "Rel_Abundance") %>%
  left_join(tss_coral$sample_table %>% rownames_to_column("SampleID") %>%
              dplyr::select(SampleID, SPP, Year, Site), by = "SampleID") %>%
  inner_join(plot_asv %>% dplyr::select(ASV_ID, Species, ASV_Order, Figure_Label),
             by = c("ASV_ID", "SPP" = "Species")) %>%
  arrange(SPP, Year, Site) %>%
  group_by(SPP) %>%
  mutate(Sample_Order = as.numeric(factor(SampleID, levels = unique(SampleID)))) %>%
  ungroup()

# =============================================================================
# SECTION 4: Panels (heatmap | threshold tiles | mean abundance bar)
# =============================================================================

heat_cols <- c("white", "#FFF5EB", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#990000")
heat_vals <- scales::rescale(c(0, 0.5, 1, 2, 5, 10, 20, 40))
thr_cols  <- c("80%" = "#A6CEE3", "90%" = "#1F78B4", "100%" = "#08306B")

build_species_row <- function(sp, is_last) {
  info <- plot_asv %>% filter(Species == sp)
  lvl  <- rev(info$Figure_Label[order(info$ASV_Order)])
  hd   <- heat_long %>% filter(SPP == sp) %>% mutate(Figure_Label = factor(Figure_Label, levels = lvl))
  info <- info %>% mutate(Figure_Label = factor(Figure_Label, levels = lvl))

  p_heat <- ggplot(hd, aes(x = reorder(SampleID, Sample_Order), y = Figure_Label, fill = Rel_Abundance)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat_cols, values = heat_vals, limits = c(0, NA)) +
    labs(x = NULL, y = NULL, title = sp) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14, face = "bold.italic"),
          legend.position = "none", plot.margin = margin(2, 0, 2, 2))

  thr <- info %>%
    dplyr::select(Figure_Label, Core_80, Core_90, Core_100) %>%
    pivot_longer(-Figure_Label, names_to = "Threshold", values_to = "Present") %>%
    mutate(Threshold = factor(sub("Core_", "", Threshold) %>% paste0("%"), levels = names(thr_cols))) %>%
    filter(Present)

  p_thr <- ggplot(thr, aes(x = Threshold, y = Figure_Label, color = Threshold)) +
    geom_point(shape = 15, size = 4.5) +
    scale_color_manual(values = thr_cols) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE, limits = lvl) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          panel.grid = element_blank(), legend.position = "none", plot.margin = margin(2, 2, 2, 0))

  max_ab <- max(info$Mean_Relative_Abundance_Percent)
  p_ab <- ggplot(info, aes(x = Mean_Relative_Abundance_Percent, y = Figure_Label)) +
    geom_col(fill = species_colors[[sp]], alpha = 0.8, width = 0.7) +
    geom_text(aes(x = Mean_Relative_Abundance_Percent + max_ab * 0.04, label = shared_symbol),
              size = 5, fontface = "bold", hjust = 0, vjust = 0.75) +
    scale_x_continuous(limits = c(0, max_ab * 1.25), expand = c(0, 0)) +
    labs(x = if (is_last) "Mean relative abundance (%)" else NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = "bold"),
          panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 0))

  plot_grid(p_heat, p_thr, p_ab, nrow = 1, rel_widths = c(0.68, 0.12, 0.20), align = "h", axis = "tb")
}

sp_present <- target_species[target_species %in% plot_asv$Species]
rows <- lapply(seq_along(sp_present), function(i) build_species_row(sp_present[i], i == length(sp_present)))
n_rows <- plot_asv %>% count(Species) %>% arrange(match(Species, sp_present)) %>% pull(n)
p_body <- plot_grid(plotlist = rows, ncol = 1, rel_heights = n_rows + 2.5, align = "v")

# --- legends -------------------------------------------------------------------
legend_counts <- core_counts %>%
  group_by(Species) %>%
  summarise(label = paste0(dplyr::first(Species), " (", paste(n[order(Threshold)], collapse = ", "), ")"),
            .groups = "drop")
sp_labels <- setNames(legend_counts$label, legend_counts$Species)

leg_heat <- get_legend(
  ggplot(data.frame(x = 1:8, f = c(0, 0.5, 1, 2, 5, 10, 20, 40)), aes(x, x, fill = f)) + geom_tile() +
    scale_fill_gradientn(colors = heat_cols, values = heat_vals, breaks = c(0, 5, 10, 20, 40),
                         name = "Relative\nabundance (%)") +
    theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)))
leg_thr <- get_legend(
  ggplot(data.frame(x = 1:3, t = factor(names(thr_cols), levels = names(thr_cols))), aes(x, 1, color = t)) +
    geom_point(shape = 15, size = 5) + scale_color_manual(values = thr_cols, name = "Prevalence") +
    theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)))
leg_sp <- get_legend(
  ggplot(data.frame(x = 1:4, s = factor(target_species, levels = target_species)), aes(x, 1, fill = s)) +
    geom_tile() + scale_fill_manual(values = species_colors, labels = sp_labels,
                                    name = "Core ASVs (80, 90, 100%)") +
    theme(legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12, face = "italic")))

legends <- plot_grid(NULL, leg_heat, leg_thr, leg_sp, ncol = 1, rel_heights = c(0.6, 0.4, 0.25, 0.35))
p_final <- plot_grid(p_body, legends, nrow = 1, rel_widths = c(0.8, 0.2))

fig_h <- max(7, 0.28 * sum(n_rows + 2.5))
ggsave(file.path(out_dir, "Figure3_core_microbiome_1pct.png"), p_final,
       width = 14, height = fig_h, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure3_core_microbiome_1pct.pdf"), p_final,
       width = 14, height = fig_h, bg = "white")

cat("\nFigure 3 saved (", sum(n_rows), "rows,", round(fig_h, 1), "in tall).\n")
cat("Appendix 2 written:", nrow(appendix2), "ASV x species rows.\n")
