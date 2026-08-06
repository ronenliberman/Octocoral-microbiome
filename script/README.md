# Florida Octocoral Microbiome Analysis — Script Folder

## Overview

This folder now contains the revised manuscript analysis scripts for the Florida octocoral microbiome project. The scripts implement the current revision workflow for core microbiome analysis, differential abundance testing, and targeted phylogenetic investigations.

## Script Folder Contents

- `Octocoral_Microbiome_Comprehensive_Analysis.Rmd` — main analysis pipeline
- `Community compsition ASV_top 20 percent per sample_MARCH17.R`
- `Core_Microbiome_Revised_AllThresholds.R`
- `Endozoicomonas_Phylogenetic_Analysis_FINAL.R`
- `Figure4_species_ANCOMBC.R`
- `Mycoplasma_Phylogenetic_Analysis.R`
- `NCBI_BLAST_unclassified_ASVs_localBLAST.R`
- `Octocoral_Revised_Analysis_COMPLETE.R`
- `Supplementary_Materials.R`

## Notes

- The prior script folder contents have been removed as part of this revision.
- The current scripts are the updated revision versions and should be used for the current manuscript submission.
- `Octocoral_Microbiome_Comprehensive_Analysis.Rmd` is the primary pipeline document and should be run first.

## How to Run

### Option 1: Run in RStudio
1. Open `Octocoral_Microbiome_Comprehensive_Analysis.Rmd`
2. Update any local file paths at the top of the document
3. Knit the document to generate the full analysis output

### Option 2: Render from the command line

```r
Rscript -e "rmarkdown::render('script/Octocoral_Microbiome_Comprehensive_Analysis.Rmd')"
```

## Recommended Order

1. `Octocoral_Microbiome_Comprehensive_Analysis.Rmd`
2. `Core_Microbiome_Revised_AllThresholds.R`
3. `Figure4_species_ANCOMBC.R`
4. `Endozoicomonas_Phylogenetic_Analysis_FINAL.R`
5. `Mycoplasma_Phylogenetic_Analysis.R`
6. `NCBI_BLAST_unclassified_ASVs_localBLAST.R`
7. `Supplementary_Materials.R`

## Dependencies

The analysis uses the same package suite described in the repository root README, including `microeco`, `qiime2R`, `vegan`, `phyloseq`, `ANCOMBC`, and visualization packages.

## Status

The folder has been updated to the current manuscript revision files. Old files present before this update have been removed.

4. **Module 4: Differential Abundance** (Lines 642-850)
   - ANCOM-BC2 analysis
   - Host-specific associates
   - Effect size estimation

5. **Module 5: Core Microbiome** (Lines 852-1150)
   - Multi-threshold analysis
   - Bootstrap validation
   - Robust core identification

6. **Module 6: Data Export** (Lines 1152-1280)
   - Results compilation
   - File generation
   - Summary statistics

## Troubleshooting

### Common Issues

**Problem**: Package installation errors
```r
# Solution: Install dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
```

**Problem**: Memory limitations
```r
# Solution: Increase memory limit or reduce rarefaction depth
memory.limit(size = 32000)  # Windows only
# Or reduce sample.size in rarefaction
```

**Problem**: File path errors
```r
# Solution: Use absolute paths or check working directory
setwd("~/path/to/Octocoral-microbiome")
list.files("raw data/")  # Verify files are present
```

**Problem**: QIIME2 import fails
```r
# Solution: Verify QIIME2 artifacts are valid
# Reimport from QIIME2 if needed
qiime tools peek table-NoEMC.qza  # In terminal
```

### Getting Help
- Check detailed documentation: `Octocoral_Analysis_Documentation.md`
- Review microeco package documentation: https://chiliubio.github.io/microeco_tutorial/
- Open an issue on GitHub: [Issues page](https://github.com/ronenliberman/Octocoral-microbiome/issues)

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Laboratory**: Joe Lopez Lab, Halmos College of Arts and Sciences, Nova Southeastern University
- **Collaborators**: 
  - Matias Gomez-Corales (MGC) - Phylogenetic analysis
  - Jose Victor Lopez (JVL) - Supervision and conceptualization
- **Sequencing**: Lopez Lab (NSU)
- **Computing Resources**: NSU High Performance Computing
- **Field Support**: [Field assistance acknowledgments]

## Contact

**Ronen Liberman**
- Email: rliberma@nova.edu
- GitHub: [@ronenliberman](https://github.com/ronenliberman)
- Lab: [Lopez Lab, NSU](https://www.nova.edu/ocean/lopez/)


### Methodological References
- Callahan et al. (2016) - DADA2 for high-resolution sample inference
- Lin & Peddada (2020) - ANCOM-BC2 methodology
- Liu et al. (2021) - microeco R package
- Bolyen et al. (2019) - QIIME2 ecosystem

## Project Status

🟢 **Active Development**
- Main analysis pipeline: Complete ✅
- Bootstrap validation: Complete ✅
- Species-specific ordinations: Complete ✅
- Manuscript: In preparation 📝
- Peer review: Pending submission

---

**Keywords**: octocoral, microbiome, 16S rRNA, coral reef, Florida, bacterial communities, core microbiome, symbiosis, holobiont, ANCOM-BC2, bootstrap validation

**Last Updated**: December 2025
