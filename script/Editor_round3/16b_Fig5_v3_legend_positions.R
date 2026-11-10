# =============================================================================
# 16b_Fig5_v3_legend_positions.R — Fig 5 v3b (chosen 20 Oct 2026)
# Same as 16_tree_figures_v3.R Figure 5, except the legend box is lower: anchored
# by its bottom-right corner, ~1 in above the Kistimonas name.
# Writes Figure5_..._tree_v3b.png; Figure 6 v3 (16_tree_figures_v3.R) is final.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

fig5_h <- 16
fig5_legend_y <- 0.175                # box bottom edge (npc); Kistimonas label top ~ 0.12
fig5_layout <- list(xlim_mult = 1.17, outgroup_above = TRUE, legend_scale = 1.5, legend_text_scale = 1.275,
                    legend_box = TRUE, host_nrow = 2, source_nrow = 1, key_gap_cm = 0.6,
                    legend_x = 0.95, legend_just = c(1, 0))

e <- tree_tip_data("Endo_r3", "Endo_all_ASVs_with_status.csv")
p <- plot_ml_tree_support_values(e$tree, e$tip_data, outgroup_pattern = "Kistimonas",
                                 fig_w = 11, fig_h = fig5_h, legend_y = fig5_legend_y,
                                 tag = "Figure 5 v3b", layout = fig5_layout)
f <- file.path(out_dir, "Figure5_Endozoicomonadaceae_tree_v3b.png")
ggsave(f, p, width = 11, height = fig5_h, dpi = 600, bg = "white"); cat("saved:", f, "\n")
