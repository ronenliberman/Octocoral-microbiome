################################################################################
# Octocoral Microbiome Revised Analysis - COMPLETE CONSOLIDATED SCRIPT
# Author: Ronen Liberman
# Date: March 2026
# Description: Full analysis pipeline for manuscript revision
#   Addresses reviewer comments on normalization, ASV-level consistency,
#   NCBI validation, CLR/Aitchison distance, species-specific ANCOM-BC2
#
# DATA FLOW:
#   d (unrarefied) --> mt_rarefied (alpha diversity)
#                  --> mt_TSS (community composition, relative abundance)
#                  --> mt_CLR (beta diversity, Aitchison distance)
#   ANCOM-BC2 uses d (unrarefied counts)
################################################################################

# =============================================================================
# SECTION 0: PACKAGES AND SETUP
# =============================================================================

library(tidyverse)
library(microeco)
library(file2meco)
library(qiime2R)
library(Biostrings)
library(vegan)
library(ggplot2)
library(cowplot)
library(stringr)
library(ANCOMBC)
library(dunn.test)
library(ggsignif)
library(flextable)
library(officer)
library(RColorBrewer)
library(magrittr)
library(scales)

setwd(here::here())   # repository root (found by the .Rproj file)
output_dir <- file.path(here::here(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

# --- Consistent colors and theme ---
species_colors <- c(
  "B. asbestinum" = "#E41A1C", "E. caribaeorum" = "#FF7F00",
  "E. flexuosa" = "#984EA3", "M. muricata" = "#4DAF4A", "Seawater" = "#377EB8"
)

my_theme <- theme_classic() + theme(
  axis.title = element_text(size = 16, face = "bold"),
  axis.text = element_text(size = 14, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  strip.text = element_text(size = 14, face = "bold.italic"),
  strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
  legend.position = "right"
)

species_list <- c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")

# =============================================================================
# SECTION 1: LOAD DATA AND FILTER
# =============================================================================

cat("\n=== SECTION 1: Loading and filtering data ===\n")

data <- qiime2meco(
  "raw files/table-NoEMC-PEE2.qza",
  taxonomy_table = "raw files/taxonomy-PEE2.qza",
  sample_table = "raw files/metadata_all_edited.csv",
  auto_tidy = TRUE
)

d <- clone(data)
d$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
d$sample_table %<>% subset(SPP %in% c("BAST", "ECAR", "MMUR", "EFLE", "SW"))
d$sample_table$Site <- gsub("DC2", "DC3", d$sample_table$Site)
d$sample_table %<>% subset(Site %in% c("DC1", "DC3", "BC4", "BC6"))

# Remove 3 unreliable samples
samples_to_remove <- c("ECAR_DC3_3_2_2023", "BAST_DC3_1_2024", "EFLE_DC3_3_2024")
d$sample_table <- d$sample_table %>% filter(!rownames(.) %in% samples_to_remove)

# Rename species
species.names <- list("SW"="Seawater","BAST"="B. asbestinum","ECAR"="E. caribaeorum",
                       "MMUR"="M. muricata","EFLE"="E. flexuosa")
d$sample_table <- d$sample_table %>% mutate(SPP = recode(SPP, !!!species.names))
d$tidy_dataset()
cat("Final dataset: Samples:", nrow(d$sample_table), "| ASVs:", nrow(d$otu_table), "\n")

# --- Sequencing summary for Results text ---
cat("\n--- Sequencing Summary (for manuscript) ---\n")

total_reads <- sum(d$otu_table)
total_asvs  <- nrow(d$otu_table)
total_samples <- nrow(d$sample_table)

cat("Total reads (after filtering):", format(total_reads, big.mark = ","), "\n")
cat("Total ASVs:", format(total_asvs, big.mark = ","), "\n")
cat("Total samples:", total_samples, "\n")

# Per-sample read depth
sample_depths <- colSums(d$otu_table)
cat("Reads per sample: min =", min(sample_depths), 
    "| max =", format(max(sample_depths), big.mark = ","),
    "| mean =", round(mean(sample_depths)), 
    "| median =", round(median(sample_depths)), "\n")

# Samples per species (including seawater)
cat("\nSamples per species:\n")
print(table(d$sample_table$SPP))

# Samples removed below rarefaction threshold
below_thresh <- sum(sample_depths < 20000)
cat("\nSamples below 20,000 reads:", below_thresh, "\n")
if (below_thresh > 0) {
  cat("  IDs:", paste(names(sample_depths[sample_depths < 20000]), collapse = ", "), "\n")
}
# =============================================================================
# SECTION 2: NORMALIZATION
# =============================================================================

cat("\n=== SECTION 2: Normalization ===\n")

# Rarefaction (for alpha diversity)
mt_rarefied <- trans_norm$new(dataset = d)$norm(method = "rarefy", sample.size = 20000)

# TSS (for community composition and relative abundance)
mt_TSS <- trans_norm$new(dataset = d)$norm(method = "TSS")

# CLR (for Aitchison beta diversity)
mt_CLR <- trans_norm$new(dataset = d)$norm(method = "rclr")

# =============================================================================
# SECTION 3: ALPHA DIVERSITY (SUPPLEMENTARY)
# =============================================================================
# See: Alpha diversity code block from conversation
# Uses mt_rarefied, calculates Observed + Shannon + Pielou
# Kruskal-Wallis + Dunn's post-hoc + ggsignif boxplots
# Outputs: p_richness, p_evenness, alpha_summary CSV

cat("\n=== SECTION 3: Alpha Diversity (for Supplementary) ===\n")

mt_rarefied$cal_alphadiv(PD = FALSE)

mt_rarefied$alpha_diversity$Pielou <- mt_rarefied$alpha_diversity$Shannon /
  log(mt_rarefied$alpha_diversity$Observed)

alpha_df <- mt_rarefied$alpha_diversity %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(
    mt_rarefied$sample_table %>% tibble::rownames_to_column("sample_id"),
    by = "sample_id"
  ) %>%
  filter(SPP != "Seawater")

alpha_summary <- alpha_df %>%
  group_by(SPP) %>%
  summarise(
    n = n(),
    Observed_mean = round(mean(Observed, na.rm = TRUE), 0),
    Observed_sd = round(sd(Observed, na.rm = TRUE), 0),
    Shannon_mean = round(mean(Shannon, na.rm = TRUE), 2),
    Shannon_sd = round(sd(Shannon, na.rm = TRUE), 2),
    Pielou_mean = round(mean(Pielou, na.rm = TRUE), 3),
    Pielou_sd = round(sd(Pielou, na.rm = TRUE), 3),
    .groups = "drop"
  )
print(alpha_summary)
write.csv(alpha_summary, file.path(output_dir, "alpha_diversity_summary_supplementary.csv"), row.names = FALSE)

kw_obs <- kruskal.test(Observed ~ SPP, data = alpha_df)
print(kw_obs)
dunn_obs <- dunn.test(alpha_df$Observed, alpha_df$SPP, method = "bonferroni")

kw_pielou <- kruskal.test(Pielou ~ SPP, data = alpha_df)
print(kw_pielou)
dunn_pielou <- dunn.test(alpha_df$Pielou, alpha_df$SPP, method = "bonferroni")

get_sig_brackets <- function(dunn_result, y_max) {
  data.frame(comparison = dunn_result$comparisons, p_adj = dunn_result$P.adjusted, stringsAsFactors = FALSE) %>%
    mutate(label = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "")) %>%
    filter(label != "") %>%
    tidyr::separate(comparison, into = c("group1", "group2"), sep = " - ", remove = FALSE) %>%
    mutate(y_position = y_max * (1.05 + 0.07 * (row_number() - 1)))
}
species_list <- c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")

alpha_df$SPP <- factor(alpha_df$SPP, levels = species_list)

sig_obs <- get_sig_brackets(dunn_obs, max(alpha_df$Observed, na.rm = TRUE))
p_richness <- ggplot(alpha_df, aes(x = SPP, y = Observed)) +
  geom_boxplot(fill = "grey85", outlier.shape = 21, outlier.size = 2, outlier.fill = "grey50", width = 0.6) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, face = "italic", color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.8), plot.margin = margin(20, 20, 10, 20)) +
  labs(y = "Observed Richness")
for (i in seq_len(nrow(sig_obs))) {
  p_richness <- p_richness + geom_signif(comparisons = list(c(sig_obs$group1[i], sig_obs$group2[i])),
                                         annotations = sig_obs$label[i], y_position = sig_obs$y_position[i], tip_length = 0.02, textsize = 5, vjust = -0.2)
}
print(p_richness)
ggsave(file.path(output_dir, "Alpha_diversity_Observed_Richness.png"), p_richness, dpi = 300, width = 8, height = 7, bg = "white")

sig_pielou <- get_sig_brackets(dunn_pielou, max(alpha_df$Pielou, na.rm = TRUE))
p_evenness <- ggplot(alpha_df, aes(x = SPP, y = Pielou)) +
  geom_boxplot(fill = "grey85", outlier.shape = 21, outlier.size = 2, outlier.fill = "grey50", width = 0.6) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, face = "italic", color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.8), plot.margin = margin(20, 20, 10, 20)) +
  labs(y = "Pielou's Evenness")
for (i in seq_len(nrow(sig_pielou))) {
  p_evenness <- p_evenness + geom_signif(comparisons = list(c(sig_pielou$group1[i], sig_pielou$group2[i])),
                                         annotations = sig_pielou$label[i], y_position = sig_pielou$y_position[i], tip_length = 0.02, textsize = 5, vjust = -0.2)
}
print(p_evenness)
ggsave(file.path(output_dir, "Alpha_diversity_Pielou_Evenness.png"), p_evenness, dpi = 300, width = 8, height = 7, bg = "white")

cat("Section 3 complete.\n")

# =============================================================================
# SECTION 4: BETA DIVERSITY
# =============================================================================

cat("\n=== SECTION 4: Beta Diversity ===\n")

# --- 4A: CLR Aitchison PCoA (with seawater) ---
cat("\n--- 4A: CLR Aitchison PCoA with seawater ---\n")
octo_clr_sw <- clone(mt_CLR)
octo_clr_sw$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa", "Seawater"))
octo_clr_sw$sample_table$Year <- as.factor(octo_clr_sw$sample_table$Year)
octo_clr_sw$tidy_dataset()

euc_dist_sw <- vegan::vegdist(t(octo_clr_sw$otu_table), method = "euclidean")
octo_clr_sw$beta_diversity <- list(euclidean = as.matrix(euc_dist_sw))

t_beta_clr_sw <- trans_beta$new(dataset = octo_clr_sw, group = "SPP", measure = "euclidean")
t_beta_clr_sw$cal_ordination(method = "PCoA")

p_pcoa_clr_sw <- t_beta_clr_sw$plot_ordination(
  plot_color = "SPP", plot_shape = "Year",
  plot_type = c("point", "centroid"),
  centroid_segment_linetype = 3, centroid_segment_size = 0.5,
  point_size = 4, point_alpha = 0.8
) +
  scale_color_manual(name = "Species", values = species_colors,
                     labels = c("B. asbestinum" = expression(italic("B. asbestinum")),
                                "E. caribaeorum" = expression(italic("E. caribaeorum")),
                                "E. flexuosa" = expression(italic("E. flexuosa")),
                                "M. muricata" = expression(italic("M. muricata")),
                                "Seawater" = "Seawater")) +
  scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17, "2025" = 15)) +
  my_theme + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = NULL)

print(p_pcoa_clr_sw)
ggsave(file.path(output_dir, "PCoA_Aitchison_CLR_with_seawater.png"),
       p_pcoa_clr_sw, dpi = 300, width = 10, height = 8, bg = "white")

# --- 4B: CLR Aitchison PCoA (octocoral only) ---
cat("\n--- 4B: CLR Aitchison PCoA octocoral only ---\n")
octo_clr <- clone(mt_CLR)
octo_clr$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa"))
octo_clr$sample_table$Year <- as.factor(octo_clr$sample_table$Year)
octo_clr$tidy_dataset()

euc_dist <- vegan::vegdist(t(octo_clr$otu_table), method = "euclidean")
octo_clr$beta_diversity <- list(euclidean = as.matrix(euc_dist))

t_beta_clr <- trans_beta$new(dataset = octo_clr, group = "SPP", measure = "euclidean")
t_beta_clr$cal_ordination(method = "PCoA")

p_pcoa_clr <- t_beta_clr$plot_ordination(
  plot_color = "SPP", plot_shape = "Year",
  plot_type = c("point", "centroid"),
  centroid_segment_linetype = 3, centroid_segment_size = 0.5,
  point_size = 4, point_alpha = 0.8
) +
  scale_color_manual(name = "Species", values = species_colors,
                     labels = c("B. asbestinum" = expression(italic("B. asbestinum")),
                                "E. caribaeorum" = expression(italic("E. caribaeorum")),
                                "E. flexuosa" = expression(italic("E. flexuosa")),
                                "M. muricata" = expression(italic("M. muricata")))) +
  scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17, "2025" = 15)) +
  my_theme + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = NULL)

print(p_pcoa_clr)
ggsave(file.path(output_dir, "PCoA_Aitchison_CLR.png"),
       p_pcoa_clr, dpi = 300, width = 10, height = 8, bg = "white")

# --- 4C: Rarefied Bray-Curtis PCoA ---
cat("\n--- 4C: Rarefied Bray-Curtis PCoA ---\n")
octo_rare <- clone(mt_rarefied)
octo_rare$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa"))
octo_rare$sample_table$Year <- as.factor(octo_rare$sample_table$Year)
octo_rare$tidy_dataset()
octo_rare$cal_betadiv(unifrac = FALSE)

t_beta_spp <- trans_beta$new(dataset = octo_rare, group = "SPP", measure = "bray")
t_beta_spp$cal_ordination(method = "PCoA")

p_pcoa_rare <- t_beta_spp$plot_ordination(
  plot_color = "SPP", plot_shape = "Year",
  plot_type = c("point", "centroid"),
  centroid_segment_linetype = 3, centroid_segment_size = 0.5,
  point_size = 4, point_alpha = 0.8
) +
  scale_color_manual(name = "Species", values = species_colors,
                     labels = c("B. asbestinum" = expression(italic("B. asbestinum")),
                                "E. caribaeorum" = expression(italic("E. caribaeorum")),
                                "E. flexuosa" = expression(italic("E. flexuosa")),
                                "M. muricata" = expression(italic("M. muricata")))) +
  scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17, "2025" = 15)) +
  my_theme + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = NULL)

print(p_pcoa_rare)
ggsave(file.path(output_dir, "PCoA_BrayCurtis_rarefied.png"),
       p_pcoa_rare, dpi = 300, width = 10, height = 8, bg = "white")

# --- 4D: TSS Bray-Curtis PCoA ---
cat("\n--- 4D: TSS Bray-Curtis PCoA ---\n")
octo_tss <- clone(mt_TSS)
octo_tss$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa"))
octo_tss$sample_table$Year <- as.factor(octo_tss$sample_table$Year)
octo_tss$tidy_dataset()
octo_tss$cal_betadiv(unifrac = FALSE)

t_beta_tss <- trans_beta$new(dataset = octo_tss, group = "SPP", measure = "bray")
t_beta_tss$cal_ordination(method = "PCoA")

p_pcoa_tss <- t_beta_tss$plot_ordination(
  plot_color = "SPP", plot_shape = "Year",
  plot_type = c("point", "centroid"),
  centroid_segment_linetype = 3, centroid_segment_size = 0.5,
  point_size = 4, point_alpha = 0.8
) +
  scale_color_manual(name = "Species", values = species_colors,
                     labels = c("B. asbestinum" = expression(italic("B. asbestinum")),
                                "E. caribaeorum" = expression(italic("E. caribaeorum")),
                                "E. flexuosa" = expression(italic("E. flexuosa")),
                                "M. muricata" = expression(italic("M. muricata")))) +
  scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17, "2025" = 15)) +
  my_theme + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = NULL)

print(p_pcoa_tss)
ggsave(file.path(output_dir, "PCoA_BrayCurtis_TSS.png"),
       p_pcoa_tss, dpi = 300, width = 10, height = 8, bg = "white")

# --- 4E: NMDS Aitchison with seawater ---
cat("\n--- 4E: NMDS Aitchison with seawater ---\n")
octo_clr_sw2 <- clone(mt_CLR)
octo_clr_sw2$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa", "Seawater"))
octo_clr_sw2$sample_table$Year <- as.factor(octo_clr_sw2$sample_table$Year)
octo_clr_sw2$tidy_dataset()

euc_dist_nmds <- vegdist(t(octo_clr_sw2$otu_table), method = "euclidean")
set.seed(123)
nmds_result <- metaMDS(euc_dist_nmds, k = 2, trymax = 200, autotransform = FALSE)
cat("NMDS stress:", round(nmds_result$stress, 4), "\n")

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)
meta_nmds <- octo_clr_sw2$sample_table %>% tibble::rownames_to_column("Sample")
nmds_df <- nmds_scores %>% left_join(meta_nmds, by = "Sample")

centroids <- nmds_df %>%
  group_by(SPP) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

p_nmds_sw <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = SPP, shape = Year)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(data = centroids, aes(x = NMDS1, y = NMDS2, color = SPP),
             shape = 3, size = 6, stroke = 2, inherit.aes = FALSE) +
  scale_color_manual(name = "Species", values = species_colors,
                     labels = c("B. asbestinum" = expression(italic("B. asbestinum")),
                                "E. caribaeorum" = expression(italic("E. caribaeorum")),
                                "E. flexuosa" = expression(italic("E. flexuosa")),
                                "M. muricata" = expression(italic("M. muricata")),
                                "Seawater" = "Seawater")) +
  scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17, "2025" = 15)) +
  annotate("text", x = Inf, y = -Inf,
           label = paste("Stress =", round(nmds_result$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 4.5, fontface = "italic") +
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right", plot.margin = margin(10, 10, 10, 10)) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  labs(title = NULL)

print(p_nmds_sw)
ggsave(file.path(output_dir, "NMDS_Aitchison_CLR_with_seawater.png"),
       p_nmds_sw, dpi = 300, width = 10, height = 8, bg = "white")

# --- 4F: NMDS Aitchison octocoral only ---
cat("\n--- 4F: NMDS Aitchison octocoral only ---\n")
octo_clr_noSW <- clone(mt_CLR)
octo_clr_noSW$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa"))
octo_clr_noSW$sample_table$Year <- as.factor(octo_clr_noSW$sample_table$Year)
octo_clr_noSW$tidy_dataset()

euc_dist_noSW <- vegdist(t(octo_clr_noSW$otu_table), method = "euclidean")
set.seed(123)
nmds_result_noSW <- metaMDS(euc_dist_noSW, k = 2, trymax = 200, autotransform = FALSE)
cat("NMDS stress (no SW):", round(nmds_result_noSW$stress, 4), "\n")

nmds_scores2 <- as.data.frame(scores(nmds_result_noSW, display = "sites"))
nmds_scores2$Sample <- rownames(nmds_scores2)
meta_nmds2 <- octo_clr_noSW$sample_table %>% tibble::rownames_to_column("Sample")
nmds_df2 <- nmds_scores2 %>% left_join(meta_nmds2, by = "Sample")

centroids2 <- nmds_df2 %>%
  group_by(SPP) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

p_nmds <- ggplot(nmds_df2, aes(x = NMDS1, y = NMDS2, color = SPP, shape = Year)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(data = centroids2, aes(x = NMDS1, y = NMDS2, color = SPP),
             shape = 3, size = 6, stroke = 2, inherit.aes = FALSE) +
  scale_color_manual(name = "Species", values = species_colors,
                     labels = c("B. asbestinum" = expression(italic("B. asbestinum")),
                                "E. caribaeorum" = expression(italic("E. caribaeorum")),
                                "E. flexuosa" = expression(italic("E. flexuosa")),
                                "M. muricata" = expression(italic("M. muricata")))) +
  scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17, "2025" = 15)) +
  annotate("text", x = Inf, y = -Inf,
           label = paste("Stress =", round(nmds_result_noSW$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 4.5, fontface = "italic") +
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right", plot.margin = margin(10, 10, 10, 10)) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  labs(title = NULL)

print(p_nmds)
ggsave(file.path(output_dir, "NMDS_Aitchison_CLR_global.png"),
       p_nmds, dpi = 300, width = 10, height = 8, bg = "white")

cat("Section 4 complete.\n")

# =============================================================================
# SECTION 5: PERMANOVA AND BETA DISPERSION
# =============================================================================

cat("\n=== SECTION 5: PERMANOVA ===\n")

# --- 5A: Rarefied Bray-Curtis ---
cat("\n--- 5A: Rarefied Bray-Curtis PERMANOVA ---\n")
t_beta_spp$cal_manova(manova_all = TRUE)
cat("PERMANOVA (rarefied Bray-Curtis):\n")
print(t_beta_spp$res_manova)

t_beta_spp$cal_betadisper()
cat("Beta dispersion (rarefied Bray-Curtis):\n")
print(t_beta_spp$res_betadisper)

# --- 5B: TSS Bray-Curtis ---
cat("\n--- 5B: TSS Bray-Curtis PERMANOVA ---\n")
t_beta_tss$cal_manova(manova_all = TRUE)
cat("PERMANOVA (TSS Bray-Curtis):\n")
print(t_beta_tss$res_manova)

# --- 5C: CLR Aitchison ---
cat("\n--- 5C: CLR Aitchison PERMANOVA ---\n")
t_beta_clr$cal_manova(manova_all = TRUE)
cat("PERMANOVA (CLR Aitchison):\n")
print(t_beta_clr$res_manova)

# --- 5D: Multi-factor PERMANOVA on Aitchison distance (primary analysis) ---
cat("\n--- 5D: Multi-factor PERMANOVA (Aitchison) ---\n")

ait_dist <- as.dist(octo_clr$beta_diversity$euclidean)
meta_clr <- octo_clr$sample_table

# Species × Year
permanova_ait_spp_year <- adonis2(
  ait_dist ~ SPP * Year, data = meta_clr, permutations = 999, by = "terms"
)
cat("PERMANOVA Aitchison: Species × Year\n")
print(permanova_ait_spp_year)

# Species × Site
permanova_ait_spp_site <- adonis2(
  ait_dist ~ SPP * Site, data = meta_clr, permutations = 999, by = "terms"
)
cat("PERMANOVA Aitchison: Species × Site\n")
print(permanova_ait_spp_site)

# Beta dispersion
t_beta_clr$cal_betadisper()
cat("Beta dispersion (Aitchison):\n")
print(t_beta_clr$res_betadisper)

# PCoA variance explained
eig <- t_beta_clr$res_ordination$eig
cat("\nPCoA axis variance (Aitchison):\n")
cat("  Axis 1:", round(eig[1], 1), "%\n")
cat("  Axis 2:", round(eig[2], 1), "%\n")
cat("  Axes 1+2:", round(eig[1] + eig[2], 1), "%\n")

# --- 5E: Species-specific PERMANOVA (Aitchison) ---
cat("\n--- 5E: Species-specific PERMANOVA (Aitchison) ---\n")

spp_permanova_ait <- list()

for (spp in species_list) {
  cat("\n  Processing:", spp, "\n")
  spp_clr <- clone(mt_CLR)
  spp_clr$sample_table <- subset(spp_clr$sample_table, SPP == spp)
  spp_clr$sample_table$Year <- as.factor(spp_clr$sample_table$Year)
  spp_clr$tidy_dataset()
  
  spp_ait <- vegan::vegdist(t(spp_clr$otu_table), method = "euclidean")
  spp_meta <- spp_clr$sample_table
  
  if (length(unique(spp_meta$Year)) > 1) {
    perm_yr <- adonis2(spp_ait ~ Year, data = spp_meta, permutations = 999)
    cat("    Year: R2 =", round(perm_yr$R2[1], 4), "F =", round(perm_yr$F[1], 2),
        "p =", perm_yr$`Pr(>F)`[1], "\n")
    spp_permanova_ait[[paste0(spp, "_Year")]] <- perm_yr
  }
  if (length(unique(spp_meta$Site)) > 1) {
    perm_st <- adonis2(spp_ait ~ Site, data = spp_meta, permutations = 999)
    cat("    Site: R2 =", round(perm_st$R2[1], 4), "F =", round(perm_st$F[1], 2),
        "p =", perm_st$`Pr(>F)`[1], "\n")
    spp_permanova_ait[[paste0(spp, "_Site")]] <- perm_st
  }
}

# --- 5E: Species-specific PERMANOVA (Year and Site within each species) ---
cat("\n--- 5E: Species-specific PERMANOVA ---\n")

spp_pcoa_plots <- list()
spp_permanova_results <- list()

for (spp in species_list) {
  cat("\n  Processing:", spp, "\n")
  
  spp_data <- clone(mt_rarefied)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == spp)
  spp_data$sample_table$Year <- as.factor(spp_data$sample_table$Year)
  spp_data$sample_table$Site <- factor(spp_data$sample_table$Site,
                                       levels = c("BC4", "BC6", "DC1", "DC3"),
                                       labels = c("Broward Shallow", "Broward Deep", "Miami Shallow", "Miami Deep")
  )
  spp_data$tidy_dataset()
  
  n_samples <- nrow(spp_data$sample_table)
  if (n_samples < 4) {
    cat("    Skipping (too few samples)\n")
    next
  }
  
  spp_data$cal_betadiv(method = "bray")
  spp_dist <- as.dist(spp_data$beta_diversity$bray)
  spp_meta <- spp_data$sample_table
  
  # Year PERMANOVA
  if (length(unique(spp_meta$Year)) > 1) {
    perm_year <- adonis2(spp_dist ~ Year, data = spp_meta, permutations = 999)
    cat("    PERMANOVA Year for", spp, ":\n")
    print(perm_year)
    spp_permanova_results[[paste0(spp, "_Year")]] <- perm_year
  }
  
  # Site PERMANOVA
  if (length(unique(spp_meta$Site)) > 1) {
    perm_site <- adonis2(spp_dist ~ Site, data = spp_meta, permutations = 999)
    cat("    PERMANOVA Site for", spp, ":\n")
    print(perm_site)
    spp_permanova_results[[paste0(spp, "_Site")]] <- perm_site
  }
  
  # Species-specific PCoA (for supplementary)
  t_spp <- trans_beta$new(dataset = spp_data, group = "Year", measure = "bray")
  t_spp$cal_ordination(method = "PCoA")
  
  p <- t_spp$plot_ordination(
    plot_color = "Year", plot_shape = "Site",
    plot_type = c("point", "centroid")
  ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
      legend.position = "none"
    ) +
    labs(title = spp)
  
  spp_pcoa_plots[[spp]] <- p
}

# Combine species PCoA panels
if (length(spp_pcoa_plots) > 0) {
  combined_spp_pcoa <- plot_grid(plotlist = spp_pcoa_plots, ncol = 2, align = "hv")
  ggsave(file.path(output_dir, "Species_specific_PCoA_panels.png"),
         combined_spp_pcoa, dpi = 300, width = 14, height = 12, bg = "white")
  print(combined_spp_pcoa)
}

# --- 5F: Summary table of all PERMANOVA results ---
cat("\n--- 5F: PERMANOVA Summary ---\n")

permanova_summary <- data.frame(
  Analysis = character(), Factor = character(),
  Df = numeric(), R2 = numeric(), F_stat = numeric(), P = numeric(),
  stringsAsFactors = FALSE
)

# Global results
for (label in c("Rarefied Bray-Curtis", "TSS Bray-Curtis", "CLR Aitchison")) {
  obj <- switch(label,
                "Rarefied Bray-Curtis" = t_beta_spp$res_manova,
                "TSS Bray-Curtis" = t_beta_tss$res_manova,
                "CLR Aitchison" = t_beta_clr$res_manova
  )
  permanova_summary <- rbind(permanova_summary, data.frame(
    Analysis = label, Factor = "Species",
    Df = obj$Df[1], R2 = round(obj$R2[1], 4),
    F_stat = round(obj$F[1], 2), P = obj$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  ))
}

# Multi-factor
for (i in 1:nrow(permanova_spp_year)) {
  rn <- rownames(permanova_spp_year)[i]
  if (rn %in% c("SPP", "Year", "SPP:Year")) {
    permanova_summary <- rbind(permanova_summary, data.frame(
      Analysis = "Multi-factor (SPP × Year)", Factor = rn,
      Df = permanova_spp_year[rn, "Df"],
      R2 = round(permanova_spp_year[rn, "R2"], 4),
      F_stat = round(permanova_spp_year[rn, "F"], 2),
      P = permanova_spp_year[rn, "Pr(>F)"],
      stringsAsFactors = FALSE
    ))
  }
}

for (i in 1:nrow(permanova_spp_site)) {
  rn <- rownames(permanova_spp_site)[i]
  if (rn %in% c("SPP", "Site", "SPP:Site")) {
    permanova_summary <- rbind(permanova_summary, data.frame(
      Analysis = "Multi-factor (SPP × Site)", Factor = rn,
      Df = permanova_spp_site[rn, "Df"],
      R2 = round(permanova_spp_site[rn, "R2"], 4),
      F_stat = round(permanova_spp_site[rn, "F"], 2),
      P = permanova_spp_site[rn, "Pr(>F)"],
      stringsAsFactors = FALSE
    ))
  }
}

# Species-specific
for (key in names(spp_permanova_results)) {
  res <- spp_permanova_results[[key]]
  parts <- strsplit(key, "_")[[1]]
  spp_name <- paste(parts[1:2], collapse = " ")
  factor_name <- parts[3]
  
  model_row <- which(rownames(res) == "Model")
  if (length(model_row) == 1) {
    permanova_summary <- rbind(permanova_summary, data.frame(
      Analysis = paste("Within", spp_name), Factor = factor_name,
      Df = res[model_row, "Df"],
      R2 = round(res[model_row, "R2"], 4),
      F_stat = round(res[model_row, "F"], 2),
      P = res[model_row, "Pr(>F)"],
      stringsAsFactors = FALSE
    ))
  }
}

cat("\nComplete PERMANOVA summary:\n")
print(permanova_summary)

write.csv(permanova_summary,
          file.path(output_dir, "PERMANOVA_complete_summary.csv"),
          row.names = FALSE)
cat("Saved: PERMANOVA_complete_summary.csv\n")

cat("Section 5 complete.\n")


# =============================================================================
# MAIN TEXT TABLE: PERMANOVA summary (Aitchison distance)
# =============================================================================

cat("\n--- Main Text PERMANOVA Table (Aitchison) ---\n")

main_perm_table <- data.frame(
  Factor = character(), Df = integer(), `Sum of Squares` = numeric(),
  R2 = numeric(), F_stat = numeric(), P = numeric(),
  stringsAsFactors = FALSE, check.names = FALSE
)

# SPP × Year model
for (rn in c("SPP", "Year", "SPP:Year", "Residual", "Total")) {
  if (rn %in% rownames(permanova_clr_spp_year)) {
    main_perm_table <- rbind(main_perm_table, data.frame(
      Factor = rn,
      Df = permanova_clr_spp_year[rn, "Df"],
      `Sum of Squares` = round(permanova_clr_spp_year[rn, "SumOfSqs"], 0),
      R2 = round(permanova_clr_spp_year[rn, "R2"], 4),
      F_stat = ifelse(is.na(permanova_clr_spp_year[rn, "F"]), NA,
                      round(permanova_clr_spp_year[rn, "F"], 2)),
      P = permanova_clr_spp_year[rn, "Pr(>F)"],
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
}

# Site and SPP:Site from the site model
for (rn in c("Site", "SPP:Site")) {
  if (rn %in% rownames(permanova_clr_spp_site)) {
    main_perm_table <- rbind(main_perm_table, data.frame(
      Factor = rn,
      Df = permanova_clr_spp_site[rn, "Df"],
      `Sum of Squares` = round(permanova_clr_spp_site[rn, "SumOfSqs"], 0),
      R2 = round(permanova_clr_spp_site[rn, "R2"], 4),
      F_stat = round(permanova_clr_spp_site[rn, "F"], 2),
      P = permanova_clr_spp_site[rn, "Pr(>F)"],
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
}

# Clean up factor names
main_perm_table$Factor <- recode(main_perm_table$Factor,
                                 "SPP" = "Species", "SPP:Year" = "Species × Year",
                                 "SPP:Site" = "Species × Site"
)

main_perm_table$Significance <- case_when(
  is.na(main_perm_table$P) ~ "",
  main_perm_table$P <= 0.001 ~ "***",
  main_perm_table$P <= 0.01  ~ "**",
  main_perm_table$P <= 0.05  ~ "*",
  TRUE ~ "ns"
)

main_perm_table$P <- ifelse(is.na(main_perm_table$P), "",
                            ifelse(main_perm_table$P < 0.001, "<0.001",
                                   as.character(main_perm_table$P)))

cat("\nMain text PERMANOVA table:\n")
print(main_perm_table)

write.csv(main_perm_table,
          file.path(output_dir, "Table_main_PERMANOVA_Aitchison.csv"),
          row.names = FALSE)

# --- Word table ---
ft_perm <- flextable(main_perm_table) %>%
  set_header_labels(
    Factor = "Factor", Df = "df",
    `Sum.of.Squares` = "Sum of Squares",
    R2 = "R²", F_stat = "F", P = "P", Significance = ""
  ) %>%
  theme_vanilla() %>%
  bold(~ Factor %in% c("Species", "Residual", "Total"), j = "Factor") %>%
  align(j = 2:7, align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  border_inner_h(border = fp_border(width = 0.5, color = "grey70")) %>%
  hline(i = 5, border = fp_border(width = 1.5, color = "black")) %>%
  add_footer_lines(
    paste0("PERMANOVA on Aitchison distances (Euclidean distance on CLR-transformed data), ",
           "999 permutations, sequential (Type I) sums of squares. ",
           "Species × Year and Species × Site are from separate models. ",
           "*** P < 0.001, ** P < 0.01, * P < 0.05.")
  ) %>%
  fontsize(size = 9, part = "footer") %>%
  set_caption("Table X. PERMANOVA results for factors influencing octocoral-associated bacterial community composition based on Aitchison distances.")

doc_perm <- read_docx() %>% body_add_flextable(ft_perm)
print(doc_perm, target = file.path(output_dir, "Table_main_PERMANOVA_Aitchison.docx"))

cat("Saved: Table_main_PERMANOVA_Aitchison.docx\n")

# =============================================================================
# SECTION 6: ANCOM-BC2 GLOBAL (between all 4 species)
# =============================================================================

cat("\n=== SECTION 6: ANCOM-BC2 Global ===\n")

# Prepare unrarefied data
target_global <- clone(d)
target_global$sample_table %<>%
  subset(SPP %in% species_list) %>%
  mutate(Year = as.factor(Year))
target_global$tidy_dataset()

# Rename unclassified taxa
rename_tax <- target_global$tax_table
for (i in 1:nrow(rename_tax)) {
  unclass_idx <- grep("__$", rename_tax[i, ])
  if (length(unclass_idx) > 0) {
    rename_tax[i, unclass_idx] <- paste0(
      rename_tax[i, unclass_idx],
      "unclassified ", gsub(".__", "", rename_tax[i, unclass_idx[1] - 1])
    )
  }
}
target_global$tax_table <- rename_tax

cat("Samples:", nrow(target_global$sample_table), "| ASVs:", nrow(target_global$otu_table), "\n")

ancom_global <- trans_diff$new(
  dataset = target_global, method = "ancombc2", group = "SPP",
  taxa_level = "all", filter_thres = 0.001, p_adjust_method = "fdr", alpha = 0.05
)

ancom_global_results <- ancom_global$res_diff
cat("Total ASVs tested:", length(unique(ancom_global_results$Taxa)), "\n")
cat("Significant ASVs:", sum(ancom_global_results$diff == TRUE, na.rm = TRUE), "\n")
print(table(ancom_global_results$Significance, useNA = "ifany"))

write.csv(ancom_global_results,
          file.path(output_dir, "ANCOMBC2_global_between_species_results.csv"), row.names = FALSE)

# Clean taxonomy function (used throughout sections 6-8)
clean_taxonomy <- function(taxa_names) {
  sapply(taxa_names, function(name) {
    parts <- unlist(strsplit(name, "\\|"))
    levels <- list(
      genus = parts[grepl("^g__", parts)], family = parts[grepl("^f__", parts)],
      order = parts[grepl("^o__", parts)], class = parts[grepl("^c__", parts)],
      phylum = parts[grepl("^p__", parts)]
    )
    for (level_name in names(levels)) {
      level <- levels[[level_name]]
      if (length(level) > 0 && level[1] != paste0(substr(level_name, 1, 1), "__")) {
        clean_name <- gsub("^..__", "", level[1])
        if (!grepl("unclassified", clean_name))
          return(if (level_name == "genus") clean_name else paste("unclassified", clean_name))
      }
    }
    return("unclassified Bacteria")
  }, USE.NAMES = FALSE)
}

top_global_sig <- ancom_global_results %>%
  filter(diff == TRUE) %>%
  group_by(Taxa) %>%
  summarise(min_p_adj = min(P.adj, na.rm = TRUE),
            max_abs_lfc = max(abs(lfc), na.rm = TRUE), .groups = "drop") %>%
  arrange(min_p_adj) %>% head(30) %>%
  mutate(ASV = str_extract(Taxa, "[a-f0-9]{32}$"), clean_tax = clean_taxonomy(Taxa))

cat("\nTop 30 globally significant ASVs:\n")
print(top_global_sig %>% dplyr::select(clean_tax, min_p_adj, max_abs_lfc))

# --- 6B: Cross-reference with NCBI BLAST ---
cat("\n--- 6B: NCBI BLAST cross-reference ---\n")

ncbi_blast_file <- file.path(here::here(), "blast", "NCBI_BLAST_poorly_classified_ASVs_summary.csv")
if (file.exists(ncbi_blast_file)) {
  ncbi_blast <- read.csv(ncbi_blast_file, stringsAsFactors = FALSE)
  
  top_sig_with_asv <- ancom_global_results %>%
    filter(diff == TRUE) %>%
    group_by(Taxa) %>%
    summarise(min_p_adj = min(P.adj, na.rm = TRUE),
              max_abs_lfc = max(abs(lfc), na.rm = TRUE), .groups = "drop") %>%
    arrange(min_p_adj) %>% head(30) %>%
    mutate(ASV_ID = str_extract(Taxa, "[a-f0-9]{32}$"), clean_tax = clean_taxonomy(Taxa))
  
  top_sig_ncbi <- top_sig_with_asv %>%
    left_join(
      ncbi_blast %>% dplyr::select(ASV_ID, total_reads, n_samples, best_hit_pident,
                                   silva_Kingdom, silva_Phylum, silva_Class, silva_Order, silva_Family, silva_Genus,
                                   ncbi_Kingdom, ncbi_Phylum, ncbi_Class, ncbi_Order, ncbi_Family, ncbi_Genus, ncbi_Species,
                                   starts_with("NCBI_hit1")),
      by = "ASV_ID"
    ) %>%
    mutate(
      was_poorly_classified = !is.na(best_hit_pident),
      ncbi_improved = case_when(
        is.na(best_hit_pident) ~ "Already classified by SILVA",
        !is.na(ncbi_Genus) & ncbi_Genus != "" ~ paste0("Resolved: ", ncbi_Genus, " (", ncbi_Species, ")"),
        !is.na(ncbi_Family) & ncbi_Family != "" ~ paste0("Resolved family: ", ncbi_Family),
        !is.na(ncbi_Class) & ncbi_Class != "" ~ paste0("Resolved class: ", ncbi_Class),
        TRUE ~ "Not improved"
      )
    )
  
  n_poorly <- sum(top_sig_ncbi$was_poorly_classified, na.rm = TRUE)
  n_improved <- sum(grepl("Resolved", top_sig_ncbi$ncbi_improved), na.rm = TRUE)
  cat("Poorly classified by SILVA:", n_poorly, "\n")
  cat("Improved by NCBI:", n_improved, "\n")
  
  write.csv(top_sig_ncbi,
            file.path(output_dir, "ANCOMBC2_top30_significant_ASVs_with_NCBI.csv"), row.names = FALSE)
} else {
  cat("NCBI BLAST file not found, skipping cross-reference.\n")
}

# --- 6C: Update single ASV taxonomy ---
cat("\n--- 6C: Taxonomy update ---\n")
if (exists("top_sig_ncbi")) {
  improved_asv <- top_sig_ncbi %>%
    filter(was_poorly_classified == TRUE, grepl("Resolved", ncbi_improved))
  
  if (nrow(improved_asv) > 0) {
    asv_id <- improved_asv$ASV_ID[1]
    cat("Updating", asv_id, "Class to Gammaproteobacteria (80% pident, Class level only)\n")
    d$tax_table[asv_id, "Class"] <- "c__Gammaproteobacteria"
    mt_rarefied$tax_table[asv_id, "Class"] <- "c__Gammaproteobacteria"
  }
}

cat("Section 6 complete.\n")

# =============================================================================
# SECTION 7: ANCOM-BC2 SPECIES-SPECIFIC (within each species)
# =============================================================================

cat("\n=== SECTION 7: ANCOM-BC2 Species-Specific ===\n")

species_ancombc_results <- list()

for (spp in species_list) {
  cat("\n--- Processing:", spp, "---\n")
  
  spp_data <- clone(d)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == spp)
  spp_data$sample_table$Year <- as.factor(spp_data$sample_table$Year)
  spp_data$sample_table$Site <- as.factor(spp_data$sample_table$Site)
  spp_data$tidy_dataset()
  
  n_samples <- nrow(spp_data$sample_table)
  n_years <- length(unique(spp_data$sample_table$Year))
  n_sites <- length(unique(spp_data$sample_table$Site))
  cat("  Samples:", n_samples, "| Years:", n_years, "| Sites:", n_sites, "\n")
  
  # Rename unclassified
  rename_tax_spp <- spp_data$tax_table
  for (i in 1:nrow(rename_tax_spp)) {
    unclass_idx <- grep("__$", rename_tax_spp[i, ])
    if (length(unclass_idx) > 0) {
      rename_tax_spp[i, unclass_idx] <- paste0(
        rename_tax_spp[i, unclass_idx],
        "unclassified ", gsub(".__", "", rename_tax_spp[i, unclass_idx[1] - 1])
      )
    }
  }
  spp_data$tax_table <- rename_tax_spp
  
  # Site
  if (n_sites >= 2 && n_samples >= 10) {
    cat("  ANCOM-BC2 by Site...\n")
    tryCatch({
      ancom_site <- trans_diff$new(
        dataset = spp_data, method = "ancombc2", group = "Site",
        taxa_level = "all", filter_thres = 0.001, p_adjust_method = "fdr", alpha = 0.05
      )
      site_results <- ancom_site$res_diff
      site_results$Species <- spp; site_results$Test <- "Site"
      cat("    Significant:", sum(site_results$diff == TRUE, na.rm = TRUE), "\n")
      species_ancombc_results[[paste0(spp, "_Site")]] <- site_results
    }, error = function(e) cat("    Failed:", conditionMessage(e), "\n"))
  }
  
  # Year
  if (n_years >= 2 && n_samples >= 10) {
    cat("  ANCOM-BC2 by Year...\n")
    tryCatch({
      ancom_year <- trans_diff$new(
        dataset = spp_data, method = "ancombc2", group = "Year",
        taxa_level = "all", filter_thres = 0.001, p_adjust_method = "fdr", alpha = 0.05
      )
      year_results <- ancom_year$res_diff
      year_results$Species <- spp; year_results$Test <- "Year"
      cat("    Significant:", sum(year_results$diff == TRUE, na.rm = TRUE), "\n")
      species_ancombc_results[[paste0(spp, "_Year")]] <- year_results
    }, error = function(e) cat("    Failed:", conditionMessage(e), "\n"))
  }
}

if (length(species_ancombc_results) > 0) {
  all_spp_ancombc <- bind_rows(species_ancombc_results)
  
  summary_spp <- all_spp_ancombc %>%
    group_by(Species, Test) %>%
    summarise(total_tested = n_distinct(Taxa),
              n_significant = sum(diff == TRUE, na.rm = TRUE), .groups = "drop")
  cat("\nSpecies-specific ANCOM-BC2 Summary:\n")
  print(summary_spp)
  
  write.csv(all_spp_ancombc,
            file.path(output_dir, "ANCOMBC2_species_specific_results.csv"), row.names = FALSE)
  write.csv(summary_spp,
            file.path(output_dir, "ANCOMBC2_species_specific_summary.csv"), row.names = FALSE)
}

cat("Section 7 complete.\n")

# =============================================================================
# SECTION 8: ANCOM-BC2 VISUALIZATION + TABLES
# =============================================================================

cat("\n=== SECTION 8: ANCOM-BC2 Visualization ===\n")

# --- 8A: LFC dot plots per species ---

make_lfc_plot <- function(df, spp_name, test_name, n_top = 20) {
  sig_df <- df %>%
    filter(diff == TRUE) %>%
    mutate(clean_tax = clean_taxonomy(Taxa)) %>%
    group_by(Taxa, clean_tax) %>%
    summarise(mean_lfc = mean(abs(lfc), na.rm = TRUE),
              min_padj = min(P.adj, na.rm = TRUE),
              lfc = lfc[1], .groups = "drop") %>%
    arrange(min_padj) %>% head(n_top)
  
  if (nrow(sig_df) == 0) return(NULL)
  sig_df <- sig_df %>% mutate(clean_tax = fct_reorder(clean_tax, lfc))
  
  ggplot(sig_df, aes(x = lfc, y = clean_tax)) +
    geom_point(aes(size = -log10(min_padj)), color = "#2166AC") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_size_continuous(name = expression(-log[10](p[adj])), range = c(2, 6)) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 9, face = "italic"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 13, face = "bold.italic", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          legend.position = "right") +
    labs(x = "Log fold change", y = NULL, title = spp_name,
         subtitle = paste("Effect of", test_name))
}

site_plots_lfc <- list()
year_plots_lfc <- list()

for (spp in unique(all_spp_ancombc$Species)) {
  site_df <- all_spp_ancombc %>% filter(Species == spp, Test == "Site")
  if (nrow(site_df) > 0 && sum(site_df$diff == TRUE, na.rm = TRUE) > 0) {
    site_plots_lfc[[spp]] <- make_lfc_plot(site_df, spp, "Site")
  }
  year_df <- all_spp_ancombc %>% filter(Species == spp, Test == "Year")
  if (nrow(year_df) > 0 && sum(year_df$diff == TRUE, na.rm = TRUE) > 0) {
    year_plots_lfc[[spp]] <- make_lfc_plot(year_df, spp, "Year")
  }
}

if (length(site_plots_lfc) > 0) {
  combined_site_lfc <- plot_grid(plotlist = site_plots_lfc, ncol = 2, align = "hv", labels = "AUTO")
  ggsave(file.path(output_dir, "ANCOMBC2_species_Site_lfc.png"),
         combined_site_lfc, dpi = 300, width = 16, height = 12, bg = "white")
  print(combined_site_lfc)
}

if (length(year_plots_lfc) > 0) {
  combined_year_lfc <- plot_grid(plotlist = year_plots_lfc, ncol = 2, align = "hv", labels = "AUTO")
  ggsave(file.path(output_dir, "ANCOMBC2_species_Year_lfc.png"),
         combined_year_lfc, dpi = 300, width = 16, height = 12, bg = "white")
  print(combined_year_lfc)
}

# --- 8B: Global between-species dot plot ---
cat("\n--- 8B: Global ANCOM-BC2 dot plot ---\n")

global_sig <- ancom_global_results %>%
  filter(diff == TRUE) %>%
  mutate(clean_tax = clean_taxonomy(Taxa)) %>%
  group_by(Taxa, clean_tax) %>%
  summarise(mean_abs_lfc = mean(abs(lfc), na.rm = TRUE),
            min_padj = min(P.adj, na.rm = TRUE),
            n_comparisons_sig = n(), .groups = "drop") %>%
  arrange(min_padj) %>% head(30) %>%
  mutate(clean_tax = fct_reorder(clean_tax, mean_abs_lfc))

p_global_dot <- ggplot(global_sig, aes(x = mean_abs_lfc, y = clean_tax)) +
  geom_point(aes(size = -log10(min_padj), color = n_comparisons_sig)) +
  scale_color_viridis_c(name = "# pairwise\ncomparisons\nsignificant", option = "D") +
  scale_size_continuous(name = expression(-log[10](p[adj])), range = c(2, 7)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, face = "italic"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "right") +
  labs(x = "Mean |log fold change|", y = NULL,
       title = "Top 30 differentially abundant ASVs among octocoral species")

ggsave(file.path(output_dir, "ANCOMBC2_global_dotplot.png"),
       p_global_dot, dpi = 300, width = 12, height = 10, bg = "white")
print(p_global_dot)

# --- 8C: Publication table of significant ASVs per species ---
cat("\n--- 8C: Publication table ---\n")

sig_table <- all_spp_ancombc %>%
  filter(diff == TRUE, Factors != "(Intercept)") %>%
  mutate(
    ASV_ID = str_extract(Taxa, "[a-f0-9]{32}$"),
    clean_tax = clean_taxonomy(Taxa),
    Factor_clean = case_when(
      Factors == "SiteBC6" ~ "Deerfield Deep vs Ref",
      Factors == "SiteDC1" ~ "Miami Shallow vs Ref",
      Factors == "SiteDC3" ~ "Miami Deep vs Ref",
      Factors == "Year2024" ~ "2024 vs 2023",
      Factors == "Year2025" ~ "2025 vs 2023",
      TRUE ~ Factors
    ),
    Direction = ifelse(lfc > 0, "Enriched", "Depleted"),
    lfc_round = round(lfc, 2),
    padj_round = formatC(P.adj, format = "e", digits = 2)
  ) %>%
  dplyr::select(Species, Test, Factor_clean, clean_tax, ASV_ID,
                Direction, lfc_round, padj_round) %>%
  arrange(Species, Test, Factor_clean, padj_round)

cat("Significant ASVs (excluding Intercept):", nrow(sig_table), "\n")

pub_table <- data.frame(
  `Host Species` = sig_table$Species,
  `Factor` = sig_table$Test,
  `Comparison` = sig_table$Factor_clean,
  `Taxon` = sig_table$clean_tax,
  `Direction` = sig_table$Direction,
  `Log FC` = sig_table$lfc_round,
  `Adj. P` = sig_table$padj_round,
  check.names = FALSE, stringsAsFactors = FALSE
)

ft <- flextable(pub_table) %>%
  merge_v(j = c("Host Species", "Factor")) %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  italic(j = "Taxon", part = "body") %>%
  italic(j = "Host Species", part = "body") %>%
  align(j = c("Direction", "Log FC", "Adj. P"), align = "center", part = "all") %>%
  fontsize(size = 9, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit() %>%
  width(j = "Taxon", width = 2.5) %>%
  width(j = "Comparison", width = 2) %>%
  border_inner_h(border = fp_border(width = 0.5, color = "grey70")) %>%
  set_caption("Differentially abundant ASVs within each octocoral species across spatial and temporal factors (ANCOM-BC2, FDR P < 0.05).")

print(ft)

doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = file.path(output_dir, "ANCOMBC2_species_specific_significant_ASVs_table.docx"))

write.csv(pub_table,
          file.path(output_dir, "ANCOMBC2_species_specific_significant_ASVs_table.csv"), row.names = FALSE)

cat("Section 8 complete.\n")

# =============================================================================
# SECTION 9: COMMUNITY COMPOSITION
# =============================================================================
# Use TSS-normalized data
octo_comp <- clone(mt_TSS)
octo_comp$sample_table %<>%
  subset(SPP %in% c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa"))
octo_comp$tidy_dataset()

species_list <- c("B. asbestinum", "E. caribaeorum", "M. muricata", "E. flexuosa")

for (spp in species_list) {
  
  cat("\n--- Processing:", spp, "---\n")
  
  # Subset to this species
  spp_data <- clone(octo_comp)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == spp)
  spp_data$tidy_dataset()
  
  n_samples <- nrow(spp_data$sample_table)
  n_asvs    <- nrow(spp_data$otu_table)
  cat("  Samples:", n_samples, "| ASVs:", n_asvs, "\n")
  
  # Get relative abundance matrix (already TSS, but ensure proportions)
  otu <- spp_data$otu_table
  rel_abund <- sweep(otu, 2, colSums(otu), "/") * 100
  
  # Rank ASVs by mean relative abundance across all samples
  asv_means <- sort(rowMeans(rel_abund), decreasing = TRUE)
  n_top <- 20  # show top 20 ASVs
  top_asvs <- names(asv_means)[1:min(n_top, length(asv_means))]
  
  cat("  Top", length(top_asvs), "ASVs account for",
      round(sum(asv_means[top_asvs]), 1), "% mean relative abundance\n")
  
  # Build ASV labels: "Genus_name (ASV rank)" using taxonomy
  tax <- spp_data$tax_table
  asv_labels <- sapply(top_asvs, function(asv_id) {
    row <- tax[asv_id, ]
    # Find the deepest informative classification
    genus_val  <- as.character(row[["Genus"]])
    family_val <- as.character(row[["Family"]])
    class_val  <- as.character(row[["Class"]])
    
    # Clean prefix
    clean <- function(x, prefix) {
      if (is.na(x) || x == "" || x == prefix) return(NA)
      name <- gsub("^.__", "", x)
      if (grepl("unclassified|uncultured|unknown", name, ignore.case = TRUE)) {
        return(paste0("uncl. ", gsub("unclassified |uncultured ", "", name)))
      }
      return(name)
    }
    
    label <- clean(genus_val, "g__")
    if (is.na(label)) label <- clean(family_val, "f__")
    if (is.na(label)) label <- clean(class_val, "c__")
    if (is.na(label)) label <- "unclassified"
    
    # Add short ASV hash for uniqueness
    short_hash <- substr(asv_id, 1, 6)
    paste0(label, " (", short_hash, ")")
  })
  
  # Build long-format dataframe for plotting
  top_rel <- rel_abund[top_asvs, , drop = FALSE]
  others  <- 100 - colSums(top_rel)
  
  plot_df <- top_rel %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ASV_ID") %>%
    mutate(ASV_label = asv_labels[ASV_ID]) %>%
    dplyr::select(-ASV_ID) %>%
    pivot_longer(cols = -ASV_label, names_to = "Sample", values_to = "Abundance")
  
  # Add "Others" category
  others_df <- data.frame(
    ASV_label = "Others",
    Sample    = names(others),
    Abundance = as.numeric(others),
    stringsAsFactors = FALSE
  )
  plot_df <- bind_rows(plot_df, others_df)
  
  # Order ASVs by mean abundance (Others always last)
  asv_order <- c(rev(asv_labels), "Others")
  plot_df$ASV_label <- factor(plot_df$ASV_label, levels = asv_order)
  
  # Color palette: distinct colors for top ASVs, grey for Others
  n_colors <- length(top_asvs)
  if (n_colors <= 20) {
    qual_pals <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))
    asv_colors <- qual_pals[1:n_colors]
  } else {
    asv_colors <- scales::hue_pal()(n_colors)
  }
  names(asv_colors) <- asv_labels
  asv_colors["Others"] <- "grey80"
  all_colors <- asv_colors[asv_order]
  
  # Plot
  p <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = ASV_label)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = all_colors, name = "ASV") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.title.x   = element_blank(),
      axis.text.x    = element_blank(),
      axis.ticks.x   = element_blank(),
      axis.title.y   = element_text(size = 16, face = "bold"),
      axis.text.y    = element_text(size = 14, color = "black"),
      legend.title   = element_text(size = 12, face = "bold"),
      legend.text    = element_text(size = 9),
      legend.key.size = unit(0.45, "cm"),
      legend.position = "right",
      plot.title     = element_text(size = 18, face = "bold.italic", hjust = 0.5),
      plot.margin    = margin(10, 10, 10, 10)
    ) +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
    labs(
      title = spp,
      y     = "Relative Abundance (%)"
    )
  
  print(p)
  
  spp_filename <- gsub("\\. ", "_", spp)
  #ggsave(
   # file.path(output_dir, paste0("Community_composition_ASV_", spp_filename, ".png")),
    #p, dpi = 300, width = 14, height = 7, bg = "white"
  #)
  cat("  Saved:", paste0("Community_composition_ASV_", spp_filename, ".png"), "\n")
}

cat("\n=== ASV-level community composition plots complete ===\n")

# =============================================================================
# SECTION 9B: COMBINED COMMUNITY COMPOSITION (unified legend, Year_Site facets)
# Creates: spp_plots, plot_df, genus_colors_unified (needed for Section 10)
# =============================================================================

cat("\n=== SECTION 9B: Combined composition figure ===\n")

# Dominant genera (>20% in any sample, excluding uncl. Bacteria)
top_genera_all <- c(
  "Endozoicomonas", "Gammaproteobacteria", "Rhodobacteraceae", "Oleiphilus",
  "Alphaproteobacteria", "Rhizobiaceae", "Mycobacterium",
  "Mycoplasma", "Bacilli", "Tistlia", "Defluviicoccus",
  "Methyloligellaceae", "Labrenzia", "Spiroplasma", "Blastopirellula",
  "SAR116_clade", "Candidatus_Hepatoplasma", "Cyclobacteriaceae", "Spirochaeta_2"
)

top_genera_display <- c(
  "Endozoicomonas"="Endozoicomonas", "Gammaproteobacteria"="c_Gammaproteobacteria",
  "Rhodobacteraceae"="f_Rhodobacteraceae", "Oleiphilus"="Oleiphilus",
  "Alphaproteobacteria"="c_Alphaproteobacteria", "Rhizobiaceae"="f_Rhizobiaceae",
  "Mycobacterium"="Mycobacterium", "Mycoplasma"="Mycoplasma",
  "Bacilli"="c_Bacilli", "Tistlia"="Tistlia", "Defluviicoccus"="Defluviicoccus",
  "Methyloligellaceae"="f_Methyloligellaceae", "Labrenzia"="Labrenzia",
  "Spiroplasma"="Spiroplasma", "Blastopirellula"="Blastopirellula",
  "SAR116_clade"="f_SAR116", "Candidatus_Hepatoplasma"="Candidatus Hepatoplasma",
  "Cyclobacteriaceae"="f_Cyclobacteriaceae", "Spirochaeta_2"="Spirochaeta"
)

genus_colors_unified <- c(
  "Endozoicomonas"="#1B9E77", "c_Gammaproteobacteria"="#7570B3",
  "f_Rhodobacteraceae"="#E7298A", "Oleiphilus"="#D95F02",
  "c_Alphaproteobacteria"="#377EB8", "f_Rhizobiaceae"="#66A61E",
  "Mycobacterium"="#B2DF8A", "Mycoplasma"="#E41A1C",
  "c_Bacilli"="#CAB2D6", "Tistlia"="#6A3D9A", "Defluviicoccus"="#FF7F00",
  "f_Methyloligellaceae"="#FDBF6F", "Labrenzia"="#B15928",
  "Spiroplasma"="#984EA3", "Blastopirellula"="#A6CEE3",
  "f_SAR116"="#33A02C", "Candidatus Hepatoplasma"="#A65628",
  "f_Cyclobacteriaceae"="#999999", "Spirochaeta"="#F781BF",
  "Non-dominant ASVs"="grey80"
)

# Helper function
build_species_composition <- function(mt_TSS_obj, species_name, site_levels, site_labels) {
  spp_data <- clone(mt_TSS_obj)
  spp_data$sample_table <- subset(spp_data$sample_table, SPP == species_name)
  spp_data$sample_table$Year <- as.factor(spp_data$sample_table$Year)
  spp_data$sample_table$Site <- factor(spp_data$sample_table$Site, levels = site_levels, labels = site_labels)
  spp_data$tidy_dataset()
  
  otu <- spp_data$otu_table
  tax <- spp_data$tax_table
  rel_abund <- sweep(otu, 2, colSums(otu), "/") * 100
  
  asv_genus <- sapply(rownames(tax), function(asv_id) {
    row <- tax[asv_id, ]
    clean <- function(x, prefix) {
      if (is.na(x) || x == "" || x == prefix) return(NA)
      name <- gsub("^.__", "", x)
      if (grepl("^unclassified|^uncultured|^unknown", name, ignore.case = TRUE))
        return(paste0("uncl. ", gsub("unclassified |uncultured ", "", name)))
      return(name)
    }
    label <- clean(as.character(row[["Genus"]]), "g__")
    if (is.na(label)) label <- clean(as.character(row[["Family"]]), "f__")
    if (is.na(label)) label <- clean(as.character(row[["Class"]]), "c__")
    if (is.na(label)) label <- "uncl. Bacteria"
    return(label)
  })
  
  genus_abund <- data.frame(Genus = asv_genus, rel_abund, check.names = FALSE) %>%
    group_by(Genus) %>% summarise(across(everything(), sum), .groups = "drop")
  
  meta <- spp_data$sample_table %>% tibble::rownames_to_column("Sample")
  genus_abund %>%
    pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance") %>%
    left_join(meta, by = "Sample") %>%
    mutate(Species = species_name)
}

site_levels <- c("BC4", "BC6", "DC1", "DC3")
site_labels <- c("Broward Shallow", "Broward Deep", "Miami Shallow", "Miami Deep")

all_df <- bind_rows(
  build_species_composition(mt_TSS, "B. asbestinum", site_levels, site_labels),
  build_species_composition(mt_TSS, "E. caribaeorum", site_levels, site_labels),
  build_species_composition(mt_TSS, "E. flexuosa", site_levels, site_labels),
  build_species_composition(mt_TSS, "M. muricata", site_levels, site_labels)
)

all_df <- all_df %>%
  mutate(Genus_plot = ifelse(Genus %in% top_genera_all, top_genera_display[Genus], "Non-dominant ASVs"))

plot_df <- all_df %>%
  group_by(Sample, Genus_plot, Species, Site, Year) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

genus_order <- c(rev(unname(top_genera_display)), "Non-dominant ASVs")
plot_df$Genus_plot <- factor(plot_df$Genus_plot, levels = genus_order)
plot_df$Species <- factor(plot_df$Species, levels = species_list)

dominant_genus <- plot_df %>%
  group_by(Sample) %>% slice_max(Abundance, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(Sample, dominant = Genus_plot, dom_abund = Abundance)

plot_df <- plot_df %>%
  left_join(dominant_genus, by = "Sample") %>%
  arrange(Species, Year, dominant, desc(dom_abund))

sample_order <- plot_df %>%
  distinct(Sample, Species, Year, dominant, dom_abund) %>%
  arrange(Species, Year, dominant, desc(dom_abund)) %>% pull(Sample)
plot_df$Sample <- factor(plot_df$Sample, levels = sample_order)

# Year_Site faceting
site_short <- c("Broward Shallow"="BrS", "Broward Deep"="BrD",
                "Miami Shallow"="MiS", "Miami Deep"="MiD")

plot_df <- plot_df %>%
  mutate(Site_short = site_short[as.character(Site)],
         Year_Site = paste0(Year, " - ", Site_short))

year_site_order <- c(
  "2023 - BrS","2023 - BrD","2023 - MiS","2023 - MiD",
  "2024 - BrS","2024 - BrD","2024 - MiS","2024 - MiD",
  "2025 - BrS","2025 - BrD","2025 - MiS","2025 - MiD"
)
plot_df$Year_Site <- factor(plot_df$Year_Site, levels = year_site_order)

# Build individual species plots
spp_plots <- list()
for (spp in species_list) {
  spp_df <- plot_df %>% filter(Species == spp)
  spp_df$Year_Site <- droplevels(spp_df$Year_Site)
  
  spp_sample_order <- spp_df %>%
    distinct(Sample, Year_Site, dominant, dom_abund) %>%
    arrange(Year_Site, dominant, desc(dom_abund)) %>% pull(Sample)
  spp_df$Sample <- factor(spp_df$Sample, levels = spp_sample_order)
  
  spp_plots[[spp]] <- ggplot(spp_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity", width = 1) +
    facet_wrap(~ Year_Site, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = genus_colors_unified, name = NULL, drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = element_text(size = 10, color = "black"),
          strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_rect(fill = "white", color = "black"),
          legend.position = "none", panel.spacing = unit(0.15, "lines"),
          plot.margin = margin(2, 5, 2, 5))
}

cat("spp_plots, plot_df, and genus_colors_unified created.\n")
cat("Section 9B complete. Ready for Section 10.\n")
# =============================================================================
# SECTION 10: COMBINED FIGURE 2
# Panel A: community composition, Panel B: NMDS, Panel C: richness, Panel D: evenness
# =============================================================================

cat("\n=== SECTION 10: Combined Figure 2 ===\n")

library(cowplot)

# --- Panel A: Community composition (no legend) ---
labeled_comp_plots <- list()
for (spp in species_list) {
  label <- ggdraw() + draw_label(spp, fontface = "bold.italic", size = 11, hjust = 0.5)
  labeled_comp_plots[[spp]] <- plot_grid(label, spp_plots[[spp]], ncol = 1, rel_heights = c(0.08, 1))
}

panel_A <- plot_grid(
  plotlist = labeled_comp_plots,
  ncol = 1,
  align = "v"
)

# --- Panel B: PCoA Aitchison (no seawater, no legend) + legend underneath ---
panel_B_plot <- p_pcoa_clr + theme(legend.position = "none")

pcoa_legend <- get_legend(
  p_pcoa_clr +
    guides(
      color = guide_legend(override.aes = list(size = 3.5), title = "Species", nrow = 2),
      shape = guide_legend(override.aes = list(size = 3.5), title = "Year", nrow = 1)
    ) +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.1, "cm")
    )
)

panel_B <- plot_grid(
  panel_B_plot, pcoa_legend,
  ncol = 1,
  rel_heights = c(1, 0.2)
)
# --- Panel C: Observed Richness ---
panel_C <- p_richness +
  theme(
    axis.text.x = element_text(size = 12, face = "italic", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(10, 10, 5, 10)
  )

# --- Panel D: Pielou's Evenness ---
panel_D <- p_evenness +
  theme(
    axis.text.x = element_text(size = 12, face = "italic", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(10, 10, 5, 10)
  )

# --- Right column: B on top, C and D below ---
right_column <- plot_grid(
  panel_B, panel_C, panel_D,
  ncol = 1,
  labels = c("B", "C", "D"),
  label_size = 16,
  label_fontface = "bold",
  rel_heights = c(1.2, 0.8, 0.8),
  align = "v"
)

# --- Left column: A with label ---
panel_A_labeled <- plot_grid(
  panel_A,
  labels = "A",
  label_size = 16,
  label_fontface = "bold"
)

# --- Main body: left + right ---
main_body <- plot_grid(
  panel_A_labeled, right_column,
  nrow = 1,
  rel_widths = c(0.6, 0.4)
)

# --- Composition legend: full width at bottom ---
comp_legend <- get_legend(
  ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors_unified, name = NULL) +
    guides(fill = guide_legend(ncol = 7, reverse = TRUE)) +
    theme(
      legend.text = element_text(size = 10, face = "italic"),
      legend.key.size = unit(0.5, "cm"),
      legend.position = "bottom",
      legend.justification = "center"
    )
)

# --- Final assembly ---
final_figure <- plot_grid(
  main_body, comp_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

print(final_figure)

ggsave(
  file.path(output_dir, "Figure2_combined_final_revision2.png"),
  final_figure, dpi = 600, width = 22, height = 18, bg = "white"
)
cat("Saved: Figure2_final_revision2.png\n")

cat("Section 10 complete.\n")
# =============================================================================
# SECTION 10: COMBINED FIGURE 2 — v2
# Panel A: community composition, Panel B: PCoA (no SW), Panel C: richness, Panel D: evenness
# A4 page dimensions (portrait): ~8.27 x 11.69 inches
# =============================================================================

cat("\n=== SECTION 10: Combined Figure 2 v2 ===\n")

library(cowplot)

# --- Panel A: Community composition (no legend) ---
labeled_comp_plots <- list()
for (spp in species_list) {
  label <- ggdraw() + draw_label(spp, fontface = "bold.italic", size = 10, hjust = 0.5)
  labeled_comp_plots[[spp]] <- plot_grid(label, spp_plots[[spp]], ncol = 1, rel_heights = c(0.07, 1))
}

panel_A <- plot_grid(
  plotlist = labeled_comp_plots,
  ncol = 1,
  align = "v"
)

# --- Panel B: PCoA Aitchison (no seawater, no legend) + legend underneath ---
panel_B_plot <- p_pcoa_clr + theme(legend.position = "none")

pcoa_legend <- get_legend(
  p_pcoa_clr +
    guides(
      color = guide_legend(override.aes = list(size = 3.5), title = "Species", nrow = 1),
      shape = guide_legend(override.aes = list(size = 3.5), title = "Year", nrow = 1)
    ) +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.1, "cm")
    )
)

panel_B <- plot_grid(
  panel_B_plot, pcoa_legend,
  ncol = 1,
  rel_heights = c(1, 0.2)
)

# --- Panel C: Observed Richness (compact) ---
panel_C <- p_richness +
  theme(
    axis.text.x = element_text(size = 8, face = "italic", color = "black",
                               angle = 45, hjust = 1),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

# --- Panel D: Pielou's Evenness (compact) ---
panel_D <- p_evenness +
  theme(
    axis.text.x = element_text(size = 8, face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

# --- Right column: B on top, C and D side by side below ---
alpha_row <- plot_grid(
  panel_C, panel_D,
  ncol = 2,
  labels = c("C", "D"),
  label_size = 14,
  label_fontface = "bold",
  align = "h"
)

right_column <- plot_grid(
  panel_B, alpha_row,
  ncol = 1,
  labels = c("B", ""),
  label_size = 14,
  label_fontface = "bold",
  rel_heights = c(1, 0.8)
)

# --- Left column: A with label ---
panel_A_labeled <- plot_grid(
  panel_A,
  labels = "A",
  label_size = 14,
  label_fontface = "bold"
)

# --- Main body ---
main_body <- plot_grid(
  panel_A_labeled, right_column,
  nrow = 1,
  rel_widths = c(0.58, 0.42)
)

# --- Composition legend: full width at bottom ---
comp_legend <- get_legend(
  ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors_unified, name = NULL) +
    guides(fill = guide_legend(ncol = 5, reverse = TRUE)) +
    theme(
      legend.text = element_text(size = 7, face = "italic"),
      legend.key.size = unit(0.35, "cm"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.spacing.x = unit(0.15, "cm")
    )
)

# --- Final assembly: A4 portrait ---
final_figure_v2 <- plot_grid(
  main_body, comp_legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

print(final_figure_v2)

# A4 = 8.27 x 11.69 inches, with margins ~7.5 x 10.5 usable
ggsave(
  file.path(output_dir, "Figure2_combined_final_v2_revision2.png"),
  final_figure_v2, dpi = 600, width = 12.5, height = 13.5, bg = "white"
)
cat("Saved: Figure2_combined_final_v2_revision2.png\n")

# Also save as PDF for journal submission
ggsave(
  file.path(output_dir, "Figure2_combined_final_v2_revision2.pdf"),
  final_figure_v2, width = 7.5, height = 10.5, bg = "white"
)
cat("Saved: Figure2_combined_final_v2.pdf\n")

cat("Section 10 complete.\n")

# --- Export standalone legends for manual adjustment ---

# Legend A: Community composition (genus colors)
legend_A <- get_legend(
  ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors_unified, name = NULL) +
    guides(fill = guide_legend(nrow = 4, reverse = TRUE)) +
    theme(
      legend.text       = element_text(size = 16, face = "italic"),
      legend.key.size   = unit(0.7, "cm"),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.justification = "center",
      legend.spacing.x  = unit(0.3, "cm"),
      legend.spacing.y  = unit(0.2, "cm")
    )
)

legend_A_plot <- plot_grid(legend_A)

ggsave(
  file.path(output_dir, "Figure2_legend_A_horizontal.png"),
  legend_A_plot, dpi = 600, width = 15, height = 6, bg = "white"
)
cat("Saved: Figure2_legend_A_composition.png\n")

# Legend B: NMDS (species colors + year shapes)
legend_B <- get_legend(
  p_pcoa_clr +
    guides(
      color = guide_legend(
        override.aes = list(size = 5),
        title = "Species",
        order = 1
      ),
      shape = guide_legend(
        override.aes = list(size = 5),
        title = "Year",
        order = 2
      )
    ) +
    theme(
      legend.text       = element_text(size = 14),
      legend.title      = element_text(size = 15, face = "bold"),
      legend.key.size   = unit(0.8, "cm"),
      legend.position   = "right",
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.4, "cm")
    )
)

legend_B_plot <- plot_grid(legend_B)

ggsave(
  file.path(output_dir, "Figure2_legend_B_pcoa.png"),
  legend_B_plot, dpi = 300, width = 5, height = 5, bg = "white"
)

# Legend B: PCoA - Species one row, Year second row (wide horizontal layout)
legend_B_wide <- get_legend(
  p_pcoa_clr +
    guides(
      color = guide_legend(
        override.aes = list(size = 5),
        title = "Species",
        nrow = 1,
        order = 1
      ),
      shape = guide_legend(
        override.aes = list(size = 5),
        title = "Year",
        nrow = 1,
        order = 2
      )
    ) +
    theme(
      legend.text       = element_text(size = 14),
      legend.title      = element_text(size = 15, face = "bold"),
      legend.key.size   = unit(1.0, "cm"),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.3, "cm"),
      legend.margin     = margin(10, 10, 10, 10)
    )
)

legend_B_wide_plot <- plot_grid(legend_B_wide)

ggsave(
  file.path(output_dir, "Figure2_legend_B_wide.png"),
  legend_B_wide_plot, dpi = 300, width = 8, height = 3, bg = "white"
)
cat("Saved: Figure2_legend_B_wide.png\n")

# =============================================================================
# EXTRACT PANEL A LEGEND ONLY — standalone export for manual editing
# =============================================================================
panel_A_legend_only <- get_legend(
  ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors_unified, name = NULL) +
    guides(fill = guide_legend(ncol = 5, reverse = TRUE)) +
    theme(
      legend.text = element_text(size = 7, face = "italic"),
      legend.key.size = unit(0.35, "cm"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.spacing.x = unit(0.15, "cm")
    )
)

panel_A_legend_plot <- plot_grid(panel_A_legend_only)

ggsave(
  file.path(output_dir, "Figure2_PanelA_legend_only.png"),
  panel_A_legend_plot, dpi = 600, width = 10, height = 2, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure2_PanelA_legend_only.pdf"),
  panel_A_legend_plot, width = 10, height = 2, bg = "white"
)
cat("Saved: Figure2_PanelA_legend_only.png/pdf\n")

