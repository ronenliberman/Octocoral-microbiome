# NCRI Octocoral Host Phylogenies

This repository contains a fully reproducible workflow for generating mitochondrial mtMutS phylogenetic trees for Caribbean octocorals using IQ TREE and visualization in R using ggtree.

## Species Included

- Briareum asbestinum
- Eunicea caribaeorum
- Eunicea flexuosa
- Muricea muricata

## Software Requirements

### External
- IQ TREE v3  
https://iqtree.github.io

### R Packages
- ggtree
- ape
- ggplot2
- patchwork

All R dependencies are automatically installed by the R Markdown file.

## Folder Structure

mtMuts/
├── BAST/
│ ├── BAST_mtMuts_IQTREE_Output/
│ └── BAST_mtMuts_Alignments/
├── ECAR/
├── EFLE/
└── MMUR/


Each species folder must contain:
- IQ TREE `.contree` output file  
- Annotation file with columns `old_label` and `new_label`

## How to Run

1. Open `NCRI_Octocoral_Host_Phylogenies.Rmd`
2. Edit this line to match your local directory:
base_dir <- "~/Desktop/NCRI_Octocoral_ID_2025/Octocoral_Host_ID_Clean_Sequences/mtMuts"


3. Click **Knit** to generate the full HTML report and final 4 panel tree.

## Output

- Individual rooted annotated trees
- Final combined 4 panel publication quality figure
- PDF export at 600 DPI


---




