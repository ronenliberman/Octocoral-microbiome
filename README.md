# Octocoral Microbiome Project

## Project Description

This repository contains the revised manuscript analysis workflow for the Florida octocoral microbiome project. The analysis focuses on bacterial community structure, core microbiomes, differential abundance testing, and targeted phylogenetic investigations of dominant symbionts in four octocoral host species.

### Study System

- **Location**: Broward County and Deerfield Beach reef sites, southeastern Florida, USA
- **Study Period**: 2023-2025
- **Hosts**: *Briareum asbestinum*, *Erythropodium caribaeorum*, *Eunicea flexuosa*, *Muricea muricata*
- **Sample Types**: Octocoral samples and seawater controls
- **Sites**: BC4, BC6, DC1, DC3
- **Sequencing**: Illumina MiSeq 2x300bp, 16S rRNA V4 region
- **Taxonomy**: SILVA 138.2

## Revision Update

The repository has been updated with the new R scripts from the manuscript revision source folder. These updated scripts now support the current revision workflow and are reflected in this README.

## Repository Contents

```
Octocoral-microbiome/
├── .git/
├── LICENSE
├── Octocoral-microbiome.Rproj
├── Octocoral_Manuscript_Revision_Analysis_Knowledge.md
├── README.md
├── maps/
├── phylogeny/
├── raw files/
├── output/
├── script/
│   ├── Octocoral_Microbiome_Comprehensive_Analysis.Rmd
│   ├── README.md
│   └── RarefactionCurve_Jun24_Oct24run.png
├── Community compsition ASV_top 20 percent per sample_MARCH17.R
├── Core_Microbiome_Revised_AllThresholds.R
├── Endozoicomonas_Phylogenetic_Analysis_CLEAN.R
├── Endozoicomonas_Phylogenetic_Analysis_FINAL.R
├── Figure4_species_ANCOMBC.R
├── Mycoplasma_Phylogenetic_Analysis.R
├── NCBI_BLAST_unclassified_ASVs_localBLAST.R
├── Octocoral_Revised_Analysis_COMPLETE.R
└── Supplementary_Materials.R
```

## Analysis Workflow

The new revision scripts include the following core components:

- Data import, filtering, and QC
- Alpha and beta diversity statistics
- Multi-threshold core microbiome analysis
- Bootstrap validation of the 90% core
- ANCOM-BC2 differential abundance testing
- Targeted phylogenetic analysis for *Endozoicomonas* and *Mycoplasma*

## How to Run

### Option 1: Open in RStudio
1. Open `script/Octocoral_Microbiome_Comprehensive_Analysis.Rmd`
2. Update local file paths if necessary
3. Knit to generate the full analysis report

### Option 2: Render from the command line

```r
Rscript -e "rmarkdown::render('script/Octocoral_Microbiome_Comprehensive_Analysis.Rmd')"
```

## Notes

- `Octocoral_Manuscript_Revision_Analysis_Knowledge.md` contains the revision rationale and analysis decisions.
- The new root-level R scripts are the updated revision code.
- If you want the revision scripts moved into `script/` for consistency, I can do that.

---

**Principal Investigator**: Ronen Liberman  
**Laboratory**: Joe Lopez Lab, Nova Southeastern University  
**Contact**: rliberma@nova.edu  

**Status**: Manuscript revision in progress
