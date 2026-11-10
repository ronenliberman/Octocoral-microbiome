# =============================================================================
# 16_tree_figures_v3.R — Figs 5 and 6, v3 legend (J. Lopez, 20 Oct 2026)
#   Host/Source legend as a compact outlined box placed in the tree's white space
#   (Fig 5: Host on 2 rows + Source on 1; Fig 6: Host 1 row + Source 1).
# Trees exactly as in v2 (15_tree_figures_v2.R): same tips, Kistimonas above its
# branch, support values; heights back to the _FIXED sizes since no legend strip.
# Plotting only: reads the saved trees and tip tables (tree_tip_data, 00_setup.R).
# Writes Figure5_..._tree_v3.png and Figure6_..._tree_v3.png.
# =============================================================================

source(here::here("script", "Editor_round3", "00_setup.R"))

# ---- layout (edit here) -------------------------------------------------------
fig5_h <- 16
fig5_legend_y <- 0.50                 # box centre height (npc); legend_x = box right edge
fig5_layout <- list(xlim_mult = 1.17, outgroup_above = TRUE, legend_scale = 1.5, legend_text_scale = 1.275,
                    legend_box = TRUE, host_nrow = 2, source_nrow = 1, key_gap_cm = 0.6,
                    legend_x = 0.95, legend_just = c(1, 0.5))
fig6_h <- 8
fig6_legend_y <- 0.95                 # box top (npc), top right as in Joe's example
fig6_layout <- list(legend_scale = 1.5, legend_text_scale = 1.275,
                    legend_box = TRUE, host_nrow = 1, source_nrow = 1, key_gap_cm = 0.6,
                    legend_x = 0.97, legend_just = c(1, 1))

# ---- Figure 5 -------------------------------------------------------------------
e <- tree_tip_data("Endo_r3", "Endo_all_ASVs_with_status.csv")
p <- plot_ml_tree_support_values(e$tree, e$tip_data, outgroup_pattern = "Kistimonas",
                                 fig_w = 11, fig_h = fig5_h, legend_y = fig5_legend_y,
                                 tag = "Figure 5", layout = fig5_layout)
f <- file.path(out_dir, "Figure5_Endozoicomonadaceae_tree_v3.png")
ggsave(f, p, width = 11, height = fig5_h, dpi = 600, bg = "white"); cat("saved:", f, "\n")

# ---- Figure 6 -------------------------------------------------------------------
m <- tree_tip_data("Myco_r3", "Myco_all_ASVs_with_status.csv", dagger = TRUE)
p <- plot_ml_tree_support_values(m$tree, m$tip_data, outgroup_pattern = "Acholeplasma",
                                 fig_w = 11, fig_h = fig6_h, legend_y = fig6_legend_y,
                                 tag = "Figure 6", layout = fig6_layout)
f <- file.path(out_dir, "Figure6_Mycoplasmataceae_tree_v3.png")
ggsave(f, p, width = 11, height = fig6_h, dpi = 600, bg = "white"); cat("saved:", f, "\n")
