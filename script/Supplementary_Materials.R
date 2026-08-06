################################################################################
# Supplementary Materials — Export Tables, PCoA Panels, and Seawater Analysis
# Author: Ronen Liberman
# Date: March 2026
# Description: 
#   Generates all supplementary tables and figures referenced in the Results.
#   Assumes objects from Octocoral_Revised_Analysis_COMPLETE.R are in memory:
#     d, mt_rarefied, mt_TSS, mt_CLR, species_list, species_colors, output_dir
#
# Add this section AFTER Section 10 in the main pipeline, or run standalone
# after sourcing the main script.
################################################################################

library(tidyverse)
library(vegan)
library(microeco)
library(flextable)
library(officer)
library(ggplot2)
library(cowplot)
library(magrittr)

set.seed(123)

# --- Output directory --------------------------------------------------------
supp_dir <- file.path("~/Postdocing_NSU/JOE LOPEZ Lab/Manuscript/Revision/Supplementary")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

species_list <- c("B. asbestinum", "E. caribaeorum", "E. flexuosa", "M. muricata")

species_colors <- c(
  "B. asbestinum" = "#E41A1C", "E. caribaeorum" = "#FF7F00",
  "E. flexuosa" = "#984EA3", "M. muricata" = "#4DAF4A", "Seawater" = "#377EB8"
)

cat("\n###############################################################\n")
cat("SUPPLEMENTARY MATERIALS\n")
cat("###############################################################\n")

# =============================================================================
# TABLE S1: Sample counts per species × site × year
# =============================================================================

cat("\n=== Table S1: Sample breakdown ===\n")

meta_octo <- d$sample_table %>%
  filter(SPP != "Seawater") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID")

# Site labels
site_labels_map <- c(
  "BC4" = "Deerfield Shallow", "BC6" = "Deerfield Deep",
  "DC1" = "Miami Shallow",     "DC3" = "Miami Deep"
)

meta_octo$Site_Label <- recode(meta_octo$Site, !!!site_labels_map)

sample_table_s1 <- meta_octo %>%
  group_by(SPP, Site_Label, Year) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(SPP, Site_Label)

# Add species totals
species_totals <- sample_table_s1 %>%
  group_by(SPP) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  mutate(Site_Label = "Total")

sample_table_s1 <- bind_rows(sample_table_s1, species_totals) %>%
  arrange(SPP, desc(Site_Label == "Total"))

write.csv(sample_table_s1, file.path(supp_dir, "Table_S1_sample_counts.csv"),
          row.names = FALSE)

# Word table
ft_s1 <- flextable(sample_table_s1) %>%
  set_header_labels(SPP = "Host Species", Site_Label = "Site") %>%
  merge_v(j = "SPP") %>%
  bold(j = "Total", part = "body") %>%
  bold(~ Site_Label == "Total", part = "body") %>%
  theme_vanilla() %>%
  italic(j = "SPP", part = "body") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  set_caption("Table S1. Number of samples collected per host species, reef site, and year.")

doc_s1 <- read_docx() %>% body_add_flextable(ft_s1)
print(doc_s1, target = file.path(supp_dir, "Table_S1_sample_counts.docx"))
cat("Saved: Table_S1\n")


# =============================================================================
# TABLE S2: Alpha diversity summary per species (mean ± SD)
# =============================================================================

cat("\n=== Table S2: Alpha diversity summary ===\n")

# Use mt_rarefied alpha diversity
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
    Observed_sd   = round(sd(Observed, na.rm = TRUE), 0),
    Shannon_mean  = round(mean(Shannon, na.rm = TRUE), 2),
    Shannon_sd    = round(sd(Shannon, na.rm = TRUE), 2),
    Pielou_mean   = round(mean(Pielou, na.rm = TRUE), 3),
    Pielou_sd     = round(sd(Pielou, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    Observed = paste0(Observed_mean, " ± ", Observed_sd),
    Shannon  = paste0(Shannon_mean, " ± ", Shannon_sd),
    Pielou   = paste0(Pielou_mean, " ± ", Pielou_sd)
  ) %>%
  dplyr::select(SPP, n, Observed, Shannon, Pielou)

# Kruskal-Wallis results
kw_obs <- kruskal.test(Observed ~ SPP, data = alpha_df)
kw_sha <- kruskal.test(Shannon ~ SPP, data = alpha_df)
kw_pie <- kruskal.test(Pielou ~ SPP, data = alpha_df)

cat("Kruskal-Wallis Observed: chi-sq =", round(kw_obs$statistic, 2),
    "p =", formatC(kw_obs$p.value, format = "e", digits = 2), "\n")
cat("Kruskal-Wallis Shannon: chi-sq =", round(kw_sha$statistic, 2),
    "p =", formatC(kw_sha$p.value, format = "e", digits = 2), "\n")
cat("Kruskal-Wallis Pielou: chi-sq =", round(kw_pie$statistic, 2),
    "p =", formatC(kw_pie$p.value, format = "e", digits = 2), "\n")

write.csv(alpha_summary, file.path(supp_dir, "Table_S2_alpha_diversity.csv"),
          row.names = FALSE)

# Dunn's pairwise results
library(dunn.test)
dunn_obs <- dunn.test(alpha_df$Observed, alpha_df$SPP, method = "bonferroni")
dunn_pie <- dunn.test(alpha_df$Pielou, alpha_df$SPP, method = "bonferroni")

dunn_df <- data.frame(
  Comparison = dunn_obs$comparisons,
  Observed_Z = round(dunn_obs$Z, 2),
  Observed_P_adj = formatC(dunn_obs$P.adjusted, format = "e", digits = 2),
  Pielou_Z = round(dunn_pie$Z, 2),
  Pielou_P_adj = formatC(dunn_pie$P.adjusted, format = "e", digits = 2),
  stringsAsFactors = FALSE
)

write.csv(dunn_df, file.path(supp_dir, "Table_S2b_dunn_pairwise.csv"),
          row.names = FALSE)

# Word table — alpha summary
ft_s2 <- flextable(alpha_summary) %>%
  set_header_labels(SPP = "Host Species", n = "n",
                    Observed = "Observed Richness", Shannon = "Shannon",
                    Pielou = "Pielou's Evenness") %>%
  theme_vanilla() %>%
  italic(j = "SPP", part = "body") %>%
  align(j = 2:5, align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  add_footer_lines(paste0(
    "Values are mean ± SD. All metrics differed significantly among species ",
    "(Kruskal-Wallis, P < 0.001). Rarefied to 20,000 reads per sample."
  )) %>%
  set_caption("Table S2. Alpha diversity metrics across octocoral host species.")

# Word table — Dunn's pairwise
ft_s2b <- flextable(dunn_df) %>%
  set_header_labels(
    Comparison = "Comparison",
    Observed_Z = "Z", Observed_P_adj = "Adj. P",
    Pielou_Z = "Z", Pielou_P_adj = "Adj. P"
  ) %>%
  add_header_row(
    values = c("", "Observed Richness", "Pielou's Evenness"),
    colwidths = c(1, 2, 2)
  ) %>%
  theme_vanilla() %>%
  align(j = 2:5, align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  add_footer_lines("Dunn's test with Bonferroni correction.") %>%
  set_caption("Table S2b. Pairwise alpha diversity comparisons.")

doc_s2 <- read_docx() %>%
  body_add_flextable(ft_s2) %>%
  body_add_par("") %>%
  body_add_flextable(ft_s2b)
print(doc_s2, target = file.path(supp_dir, "Table_S2_alpha_diversity.docx"))
cat("Saved: Table_S2\n")


# =============================================================================
# TABLE S3: PERMANOVA summary — all models (Aitchison)
# =============================================================================

cat("\n=== Table S3: PERMANOVA summary ===\n")

# --- Global one-way ---
octo_clr_noSW <- clone(mt_CLR)
octo_clr_noSW$sample_table %<>%
  subset(SPP %in% species_list) %>%
  mutate(Year = as.factor(Year))
octo_clr_noSW$tidy_dataset()

ait_dist <- vegan::vegdist(t(octo_clr_noSW$otu_table), method = "euclidean")
meta_clr <- octo_clr_noSW$sample_table

# Multi-factor
perm_spp_year <- adonis2(ait_dist ~ SPP * Year, data = meta_clr,
                         permutations = 999, by = "terms")
perm_spp_site <- adonis2(ait_dist ~ SPP * Site, data = meta_clr,
                         permutations = 999, by = "terms")

# Build summary
build_perm_row <- function(perm_result, analysis_label) {
  rows <- list()
  for (rn in rownames(perm_result)) {
    if (rn %in% c("Residual", "Total")) next
    rows[[length(rows) + 1]] <- data.frame(
      Analysis = analysis_label,
      Factor = rn,
      Df = perm_result[rn, "Df"],
      R2 = round(perm_result[rn, "R2"], 4),
      F_stat = round(perm_result[rn, "F"], 2),
      P = perm_result[rn, "Pr(>F)"],
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

perm_summary <- bind_rows(
  build_perm_row(perm_spp_year, "SPP × Year"),
  build_perm_row(perm_spp_site, "SPP × Site")
)

# Species-specific
for (spp in species_list) {
  spp_clr <- clone(mt_CLR)
  spp_clr$sample_table <- subset(spp_clr$sample_table, SPP == spp)
  spp_clr$sample_table$Year <- as.factor(spp_clr$sample_table$Year)
  spp_clr$tidy_dataset()
  
  spp_ait <- vegan::vegdist(t(spp_clr$otu_table), method = "euclidean")
  spp_meta <- spp_clr$sample_table
  
  if (length(unique(spp_meta$Year)) > 1) {
    res <- adonis2(spp_ait ~ Year, data = spp_meta, permutations = 999)
    perm_summary <- bind_rows(perm_summary, data.frame(
      Analysis = paste("Within", spp), Factor = "Year",
      Df = res$Df[1], R2 = round(res$R2[1], 4),
      F_stat = round(res$F[1], 2), P = res$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    ))
  }
  if (length(unique(spp_meta$Site)) > 1) {
    res <- adonis2(spp_ait ~ Site, data = spp_meta, permutations = 999)
    perm_summary <- bind_rows(perm_summary, data.frame(
      Analysis = paste("Within", spp), Factor = "Site",
      Df = res$Df[1], R2 = round(res$R2[1], 4),
      F_stat = round(res$F[1], 2), P = res$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    ))
  }
}

perm_summary$Significance <- case_when(
  perm_summary$P <= 0.001 ~ "***",
  perm_summary$P <= 0.01  ~ "**",
  perm_summary$P <= 0.05  ~ "*",
  TRUE ~ "ns"
)

write.csv(perm_summary, file.path(supp_dir, "Table_S3_PERMANOVA_Aitchison.csv"),
          row.names = FALSE)

ft_s3 <- flextable(perm_summary) %>%
  set_header_labels(F_stat = "F", P = "P-value", Significance = "Sig.") %>%
  merge_v(j = "Analysis") %>%
  theme_vanilla() %>%
  italic(~ grepl("Within", Analysis), j = "Analysis", part = "body") %>%
  align(j = 3:7, align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  add_footer_lines("PERMANOVA on Aitchison distances (999 permutations). *** P < 0.001, ** P < 0.01, * P < 0.05.") %>%
  set_caption("Table S3. PERMANOVA results across global and within-species models.")

doc_s3 <- read_docx() %>% body_add_flextable(ft_s3)
print(doc_s3, target = file.path(supp_dir, "Table_S3_PERMANOVA.docx"))
cat("Saved: Table_S3\n")


# =============================================================================
# TABLE S4: Seawater vs. Octocoral pairwise PERMANOVA
# =============================================================================

cat("\n=== Table S4: Seawater vs. octocoral comparisons ===\n")

# Include seawater in CLR dataset
octo_clr_sw <- clone(mt_CLR)
octo_clr_sw$sample_table %<>%
  subset(SPP %in% c(species_list, "Seawater")) %>%
  mutate(Year = as.factor(Year))
octo_clr_sw$tidy_dataset()

ait_dist_sw <- vegan::vegdist(t(octo_clr_sw$otu_table), method = "euclidean")
meta_sw <- octo_clr_sw$sample_table

# Pairwise PERMANOVA: each species vs seawater
sw_pairwise <- data.frame()

for (spp in species_list) {
  pair_meta <- meta_sw %>% filter(SPP %in% c(spp, "Seawater"))
  pair_samples <- rownames(pair_meta)
  pair_dist <- as.dist(as.matrix(ait_dist_sw)[pair_samples, pair_samples])
  
  res <- adonis2(pair_dist ~ SPP, data = pair_meta, permutations = 999)
  
  sw_pairwise <- bind_rows(sw_pairwise, data.frame(
    Comparison = paste(spp, "vs. Seawater"),
    Df = res$Df[1],
    R2 = round(res$R2[1], 4),
    F_stat = round(res$F[1], 2),
    P = res$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  ))
}

# Also pairwise among all octocoral species
for (i in 1:(length(species_list) - 1)) {
  for (j in (i + 1):length(species_list)) {
    sp1 <- species_list[i]
    sp2 <- species_list[j]
    pair_meta <- meta_sw %>% filter(SPP %in% c(sp1, sp2))
    pair_samples <- rownames(pair_meta)
    pair_dist <- as.dist(as.matrix(ait_dist_sw)[pair_samples, pair_samples])
    
    res <- adonis2(pair_dist ~ SPP, data = pair_meta, permutations = 999)
    
    sw_pairwise <- bind_rows(sw_pairwise, data.frame(
      Comparison = paste(sp1, "vs.", sp2),
      Df = res$Df[1],
      R2 = round(res$R2[1], 4),
      F_stat = round(res$F[1], 2),
      P = res$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    ))
  }
}

# Bonferroni correction
sw_pairwise$P_adj <- p.adjust(sw_pairwise$P, method = "bonferroni")
sw_pairwise$Significance <- case_when(
  sw_pairwise$P_adj <= 0.001 ~ "***",
  sw_pairwise$P_adj <= 0.01  ~ "**",
  sw_pairwise$P_adj <= 0.05  ~ "*",
  TRUE ~ "ns"
)

# Flag seawater comparisons
sw_pairwise$Type <- ifelse(grepl("Seawater", sw_pairwise$Comparison),
                           "vs. Seawater", "Between species")

write.csv(sw_pairwise, file.path(supp_dir, "Table_S4_pairwise_PERMANOVA.csv"),
          row.names = FALSE)

ft_s4 <- flextable(sw_pairwise %>% dplyr::select(-Type)) %>%
  set_header_labels(F_stat = "F", P = "P-value", P_adj = "Adj. P",
                    Significance = "Sig.") %>%
  theme_vanilla() %>%
  align(j = 2:7, align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit() %>%
  add_footer_lines("Pairwise PERMANOVA on Aitchison distances (999 permutations). P-values Bonferroni-corrected.") %>%
  set_caption("Table S4. Pairwise PERMANOVA comparisons among octocoral species and seawater.")

doc_s4 <- read_docx() %>% body_add_flextable(ft_s4)
print(doc_s4, target = file.path(supp_dir, "Table_S4_pairwise_PERMANOVA.docx"))

cat("Saved: Table_S4\n")
cat("\nSeawater comparison summary:\n")
print(sw_pairwise %>% filter(Type == "vs. Seawater") %>%
        dplyr::select(Comparison, R2, F_stat, P_adj, Significance))


# =============================================================================
# TABLE S5: ANCOM-BC2 species-specific significant ASVs
# (Already exported by main pipeline)
# =============================================================================

cat("\n=== Table S5: ANCOM-BC2 results ===\n")

# Check if already exists
ancom_csv <- file.path(output_dir, "Figure4_ANCOMBC2_Year_significant.csv")
ancom_site_csv <- file.path(output_dir, "Figure4_ANCOMBC2_Site_significant.csv")

if (file.exists(ancom_csv) & file.exists(ancom_site_csv)) {
  cat("  ANCOM-BC2 results already exported by Figure4_species_ANCOMBC.R:\n")
  cat("    - Figure4_ANCOMBC2_Year_significant.csv\n")
  cat("    - Figure4_ANCOMBC2_Site_significant.csv\n")
  cat("    - Figure4_ANCOMBC2_combined_summary_stats.csv\n")
  cat("  Importing and formatting as Word table...\n")
  
  year_sig <- read.csv(ancom_csv)
  site_sig <- read.csv(ancom_site_csv)
  
  # Combine and format
  combined_sig <- bind_rows(
    year_sig %>% mutate(Analysis = "Temporal"),
    site_sig %>% mutate(Analysis = "Spatial")
  ) %>%
    mutate(
      Taxon = clean_tax,
      LFC = round(lfc, 2),
      SE = round(se, 2),
      Adj_P = formatC(P.adj, format = "e", digits = 2),
      Direction = ifelse(lfc > 0, "Enriched", "Depleted")
    ) %>%
    dplyr::select(Species, Analysis, Comparison, Taxon, Direction, LFC, SE, Adj_P) %>%
    arrange(Species, Analysis, Comparison, Adj_P)
  
  write.csv(combined_sig, file.path(supp_dir, "Table_S5_ANCOMBC2_significant.csv"),
            row.names = FALSE)
  
  ft_s5 <- flextable(combined_sig) %>%
    set_header_labels(Adj_P = "Adj. P") %>%
    merge_v(j = c("Species", "Analysis")) %>%
    theme_vanilla() %>%
    italic(j = c("Species", "Taxon"), part = "body") %>%
    align(j = 5:8, align = "center", part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 9, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    autofit() %>%
    width(j = "Taxon", width = 2.2) %>%
    add_footer_lines("ANCOM-BC2 on unrarefied ASV counts, FDR-adjusted P < 0.05. LFC = log fold change, SE = standard error.") %>%
    set_caption("Table S5. Differentially abundant ASVs within each octocoral species (ANCOM-BC2).")
  
  doc_s5 <- read_docx() %>% body_add_flextable(ft_s5)
  print(doc_s5, target = file.path(supp_dir, "Table_S5_ANCOMBC2.docx"))
  cat("Saved: Table_S5\n")
  
} else {
  cat("  WARNING: ANCOM-BC2 CSVs not found in Output/\n")
  cat("  Run Figure4_species_ANCOMBC.R first, then re-run this section.\n")
}


# =============================================================================
# FIGURE S: Species-specific PCoA panels (Aitchison distance)
# =============================================================================

cat("\n=== Figure S: Species-specific PCoA (Aitchison) ===\n")

site_labels_map <- c(
  "BC4" = "DfS", "BC6" = "DfD",
  "DC1" = "MiS", "DC3" = "MiD"
)

year_shapes <- c("2023" = 16, "2024" = 17, "2025" = 15)
year_colors <- c("2023" = "#1b9e77", "2024" = "#d95f02", "2025" = "#7570b3")

spp_pcoa_list <- list()

for (spp in species_list) {
  cat("  Processing:", spp, "\n")
  
  spp_clr <- clone(mt_CLR)
  spp_clr$sample_table <- subset(spp_clr$sample_table, SPP == spp)
  spp_clr$sample_table$Year <- as.factor(spp_clr$sample_table$Year)
  spp_clr$sample_table$Site_Label <- recode(spp_clr$sample_table$Site, !!!site_labels_map)
  spp_clr$tidy_dataset()
  
  spp_ait <- vegan::vegdist(t(spp_clr$otu_table), method = "euclidean")
  
  # PCoA
  pcoa_res <- cmdscale(spp_ait, k = 2, eig = TRUE)
  var_exp <- pcoa_res$eig / sum(abs(pcoa_res$eig)) * 100
  
  pcoa_df <- as.data.frame(pcoa_res$points)
  colnames(pcoa_df) <- c("PCo1", "PCo2")
  pcoa_df$Sample <- rownames(pcoa_df)
  pcoa_df <- pcoa_df %>%
    left_join(
      spp_clr$sample_table %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        dplyr::select(Sample, Year, Site_Label),
      by = "Sample"
    )
  
  # Centroids by Year
  centroids <- pcoa_df %>%
    group_by(Year) %>%
    summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2), .groups = "drop")
  
  p <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Year, shape = Site_Label)) +
    geom_point(size = 3.5, alpha = 0.8) +
    geom_point(data = centroids, aes(x = PCo1, y = PCo2, color = Year),
               shape = 3, size = 5, stroke = 2, inherit.aes = FALSE) +
    scale_color_manual(values = year_colors) +
    scale_shape_manual(name = "Site", values = c("DfS" = 16, "DfD" = 17,
                                                  "MiS" = 15, "MiD" = 18)) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "right"
    ) +
    labs(
      title = spp,
      x = paste0("PCo1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PCo2 (", round(var_exp[2], 1), "%)")
    )
  
  spp_pcoa_list[[spp]] <- p
  
  # Save individual
  spp_fname <- gsub("\\. ", "_", spp)
  ggsave(file.path(supp_dir, paste0("FigS_PCoA_Aitchison_", spp_fname, ".png")),
         p, width = 8, height = 6, dpi = 300, bg = "white")
}

# Combined 4-panel
p_combined_pcoa <- plot_grid(
  plotlist = lapply(spp_pcoa_list, function(p) p + theme(legend.position = "none")),
  ncol = 2, align = "hv", labels = "AUTO", label_size = 18
)

# Extract shared legend
shared_legend <- get_legend(
  spp_pcoa_list[[1]] +
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(0, 0, 0, 0))
)

p_final_pcoa <- plot_grid(
  p_combined_pcoa, shared_legend,
  ncol = 1, rel_heights = c(1, 0.1)
)

ggsave(file.path(supp_dir, "FigS_Species_PCoA_Aitchison_combined.png"),
       p_final_pcoa, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigS_Species_PCoA_Aitchison_combined.pdf"),
       p_final_pcoa, width = 14, height = 12)

cat("Saved: Species-specific PCoA panels\n")


# =============================================================================
# TABLE S6: Core microbiome complete list
# (Already exported by Core_Microbiome_Revised_AllThresholds.R)
# =============================================================================

cat("\n=== Table S6: Core microbiome ===\n")

core_csv <- file.path(output_dir, "Core_Microbiome_AllThresholds_Summary_TSS.csv")

if (file.exists(core_csv)) {
  cat("  Already exported: Core_Microbiome_AllThresholds_Summary_TSS.csv\n")
  cat("  Copying to Supplementary folder...\n")
  file.copy(core_csv, file.path(supp_dir, "Table_S6_core_microbiome_all_thresholds.csv"),
            overwrite = TRUE)
} else {
  cat("  WARNING: Core microbiome CSV not found. Run Core_Microbiome_Revised_AllThresholds.R first.\n")
}

# Also copy bootstrap summary
boot_csv <- file.path(output_dir, "Robust_Core_90pct_Bootstrap_Summary.csv")
if (file.exists(boot_csv)) {
  file.copy(boot_csv, file.path(supp_dir, "Table_S6b_bootstrap_validated_core.csv"),
            overwrite = TRUE)
  cat("  Copied: Bootstrap validated core summary\n")
}


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n###############################################################\n")
cat("SUPPLEMENTARY MATERIALS COMPLETE\n")
cat("###############################################################\n")
cat("\nFiles saved to:", supp_dir, "\n\n")

cat("Tables:\n")
cat("  Table_S1_sample_counts.docx/.csv\n")
cat("  Table_S2_alpha_diversity.docx/.csv\n")
cat("  Table_S2b_dunn_pairwise.csv\n")
cat("  Table_S3_PERMANOVA.docx/.csv\n")
cat("  Table_S4_pairwise_PERMANOVA.docx/.csv  (includes seawater comparisons)\n")
cat("  Table_S5_ANCOMBC2.docx/.csv\n")
cat("  Table_S6_core_microbiome_all_thresholds.csv\n")
cat("  Table_S6b_bootstrap_validated_core.csv\n")

cat("\nFigures:\n")
cat("  FigS_Species_PCoA_Aitchison_combined.png/.pdf\n")
cat("  FigS_PCoA_Aitchison_[species].png  (4 individual panels)\n")
