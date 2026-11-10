# =============================================================================
# 17_Fig2_y_axis_title.R — add "Relative abundance (%)" to Figure 2A (21 Oct 2026)
# The final Figure 2 (GIMP-assembled PNG) has no y-axis title on panel A.
# This draws the title vertically in the existing left white margin, centred on
# panel A's axis (top "100" to bottom "0"). Same canvas (3816 x 4192 px,
# 300 dpi); no other pixels move or are resampled. Also whites out the clipped
# stray glyph left of the taxa legend (cols ~101-114, rows ~3939-3961).
# =============================================================================
library(magick)

# Input: the GIMP-assembled Figure 2 (panels from Octocoral_Revised_Analysis_COMPLETE.R and the
# community-composition script). Figures are not tracked in git, so place the file here first.
in_png  <- here::here("output", "Editor_round3", "Figure2_combined_final_GIMP.png")
out_png <- here::here("output", "Editor_round3", "Figure2_combined_final_ytitle.png")

# ---- layout (edit here; pixels of the 300 dpi image) --------------------------
title_txt  <- "Relative abundance (%)"
title_size <- 42          # em in px; panel A tick labels are the same size
title_font <- "Arial"
title_wt   <- 700         # bold, like the axis titles of panels B-D
axis_top   <- 137         # row of the "100" tick, top panel
axis_bot   <- 3680        # row of the "0" tick, bottom panel
title_x    <- 10          # left edge of the rotated title (free margin = cols 1-61)

im <- image_read(path.expand(in_png))
info <- image_info(im)

# stray clipped glyph left of the legend
im <- image_composite(im, image_blank(30, 40, "white"), offset = "+95+3930")

# render the title, trim to its ink, rotate to read bottom-to-top
lab <- image_blank(1200, 120, "none") |>
  image_annotate(title_txt, font = title_font, weight = title_wt, size = title_size,
                 color = "black", gravity = "west", location = "+10+0") |>
  image_trim() |>
  image_rotate(-90)
li <- image_info(lab)
y0 <- round((axis_top + axis_bot) / 2 - li$height / 2)
im <- image_composite(im, lab, offset = sprintf("+%d+%d", title_x, y0))

stopifnot(image_info(im)$width == info$width, image_info(im)$height == info$height)
image_write(im, path.expand(out_png), format = "png", density = "300x300")
cat("title box:", li$width, "x", li$height, "px at x =", title_x, ", y =", y0, "\nsaved:", out_png, "\n")
