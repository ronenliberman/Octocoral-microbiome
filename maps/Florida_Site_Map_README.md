# Florida Octocoral Site Map

This repository contains a fully reproducible R workflow for generating a publication quality map of octocoral sampling sites in southeast Florida.

The map includes:
- Miami and Broward sampling locations
- Shallow and deep site classification
- City reference points
- Custom label positioning for visual clarity
- North arrow and cartographic styling

---

## Contents

- `Florida_map.Rmd`  
  Main reproducible analysis and figure generation script.

- `figures/`  
  Output directory containing exported figures:
  - `map_site_plot_v4.pdf`
  - `map_site_plot_v4.png`

---

## Required R Packages

The script automatically installs and loads all dependencies:

- ggplot2  
- maps  
- tidyverse  
- ggspatial  
- ggrepel  
- pacman  

---

## How to Run

1. Open `Florida_map.Rmd` in RStudio
2. Knit to HTML
3. Figures will be automatically saved into the `figures/` folder

---

## Reproducibility

- All random placement behavior is controlled with `set.seed(123)`
- All file paths are relative to the project directory
- No user specific absolute paths are used
- Data are fully defined inside the script

---

## Author

**Matías Gómez-Corrales**  
National Coral Reef Institute  
Nova Southeastern University  

---

## License

MIT License recommended for open scientific reuse.
