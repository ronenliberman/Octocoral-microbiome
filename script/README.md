# Florida Octocoral Microbiome Analysis

## Overview

This repository contains the complete analytical pipeline and reproducible code for characterizing bacterial microbiomes associated with four octocoral species from southeastern Florida reefs. The analysis integrates 16S rRNA amplicon sequencing data to examine host-specificity, temporal stability, and core microbiome composition across *Briareum asbestinum*, *Erythropodium caribaeorum*, *Eunicea flexuosa*, and *Muricea muricata*.

## Repository Contents

```
Octocoral-microbiome/
├── script/
│   └── Octocoral_Microbiome_Comprehensive_Analysis.Rmd    # Main analysis pipeline
├── raw data/
│   ├── table-NoEMC.qza                                     # QIIME2 feature table
│   ├── taxonomy.qza                                        # QIIME2 taxonomy
│   └── metadata_all_edited.csv                             # Sample metadata
├── output/
│   └── [Analysis outputs - see below]
├── README.md                                               # This file
├── LICENSE                                                 # MIT License
└── Octocoral_Analysis_Documentation.md                     # Detailed methods documentation
```

## Data Description

### Study System
- **Location**: Broward County and Deerfield Beach reef sites, southeastern Florida, USA
- **Study Period**: 2023-2025
- **Sample Types**: Field-collected octocoral fragments + seawater controls
- **Sites**: BC4, BC6, DC1, DC3
- **Sample Size**: 200+ octocoral samples across four species

### Sequencing Details
- **Platform**: Illumina MiSeq (2x300bp paired-end)
- **Target**: 16S rRNA V4 region
- **Primers**: 515F/806R
- **Processing**: QIIME2 pipeline
- **Database**: SILVA Release 138.2
- **Rarefaction Depth**: 20,000 reads per sample

## Analytical Approach

### 1. Data Processing
- Quality control and filtering
- Rarefaction to standardized depth
- Taxonomic assignment using SILVA database
- Removal of contaminants and low-quality samples

### 2. Diversity Analysis
**Alpha Diversity**
- Shannon diversity index
- Observed ASV richness
- Simpson diversity
- Chao1 and ACE estimators
- Pielou's evenness

**Beta Diversity**
- Bray-Curtis dissimilarity
- Principal Coordinates Analysis (PCoA)
- PERMANOVA testing for group effects

### 3. Compositional Analysis
- Phylum and family-level composition
- Relative abundance visualization
- Taxonomic distribution patterns

### 4. Differential Abundance Testing
- **Method**: ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction 2)
- **Rationale**: Properly handles compositional data, zero-inflation, and sampling fraction bias
- **Output**: Host-specific bacterial associates with statistical validation

### 5. Core Microbiome Identification
- Multi-threshold prevalence analysis (75%, 90%, 100%)
- **Bootstrap validation** (1000 iterations) with fixed sample sizes
- Identification of "robust core" members (≥90% bootstrap frequency)
- Family-level aggregation and visualization

### 6. Species-Specific Patterns
- Individual PCoA ordinations per species
- Temporal variation analysis (across years)
- Spatial variation analysis (across sites)
- Sample-level pattern examination

## Software Requirements

### R Version
- R ≥ 4.0.0

### Required R Packages

**Core Microbiome Analysis**
```r
install.packages("microeco", repos = BiocManager::repositories())
install.packages("file2meco")
devtools::install_github("jbisanz/qiime2R")
```

**Data Manipulation & Visualization**
```r
install.packages(c("tidyverse", "magrittr", "ggplot2", "ggrepel", 
                   "ggpubr", "ggh4x", "grafify", "scales", 
                   "RColorBrewer", "paletteer"))
```

**Statistical Analysis**
```r
install.packages(c("vegan", "phyloseq", "MMUPHin", "multcompView", 
                   "pgirmess", "dunn.test", "performance"))
```

**Data Export**
```r
install.packages(c("writexl", "grid", "reshape2"))
```

### Additional Requirements
- QIIME2 (for initial sequence processing - not required if using provided .qza files)
- At least 16GB RAM recommended
- ~45 minutes processing time for full analysis

## Installation

### 1. Clone Repository
```bash
git clone https://github.com/ronenliberman/Octocoral-microbiome.git
cd Octocoral-microbiome
```

### 2. Install R Packages
```r
# Install all required packages
source("install_packages.R")  # If provided
# OR install individually (see Software Requirements above)
```

### 3. Set Working Directory
```r
setwd("path/to/Octocoral-microbiome")
```

## Usage

### Quick Start

**Option 1: Full Analysis in RStudio**
1. Open `script/Octocoral_Microbiome_Comprehensive_Analysis.Rmd`
2. Update file paths if needed (lines 11-12)
3. Click "Knit" to render entire analysis

**Option 2: Run from Command Line**
```r
Rscript -e "rmarkdown::render('script/Octocoral_Microbiome_Comprehensive_Analysis.Rmd')"
```

**Option 3: Interactive Section-by-Section**
```r
# Open .Rmd file in RStudio
# Run code chunks individually using Ctrl+Shift+Enter (Windows) or Cmd+Shift+Enter (Mac)
```

### Customization

**Modify Rarefaction Depth**
```r
# Line ~265
mt_rarefied <- tmp$norm(method = "rarefy", sample.size = 20000)  # Change 20000
```

**Adjust Core Microbiome Thresholds**
```r
# Line ~365
prevalence_thresholds <- c(75, 90, 100)  # Add/remove thresholds
```

**Change Bootstrap Parameters**
```r
# Lines ~475-477
bootstrap_iterations <- 1000      # Increase for more robust estimates
prevalence_threshold <- 90        # Core prevalence requirement
bootstrap_threshold <- 90         # Bootstrap frequency requirement
```

## Output Files

The analysis generates multiple output files in the `output/` directory:

### Core Microbiome Results
- `Core_Microbiome_All_Thresholds_All_Species.xlsx` - Excel workbook with separate sheets for each species × threshold combination
- `Robust_Core_90pct_Bootstrap_Summary.csv` - Bootstrap-validated core members with prevalence and abundance statistics

### Diversity Analyses
- `Alpha_Diversity_All_Samples.csv` - All alpha diversity metrics per sample
- `Beta_diversity_results/` - PERMANOVA statistics and PCoA coordinates

### Differential Abundance
- `Differential_Abundance_Top20.csv` - Significant differentially abundant taxa with statistics

### Visualizations (PNG format, 300 DPI)
- `Core_Abundance_Across_Thresholds.png` - Core microbiome trends across prevalence thresholds
- `Core_90pct_Family_18families_Custom_Colors.png` - Family-level core composition
- `PCoA_[Species]_all_labeled.png` - Species-specific ordinations (4 files)
- `significant_asvs_boxplot.png` - Differential abundance visualization
- `Alpha_Diversity_Shannon.png` - Shannon diversity by species
- `Sequencing_Depth_Distribution.png` - Read count distribution

## Reproducibility

### Ensuring Reproducibility
1. **Seed Setting**: All random processes use `set.seed(123)`
2. **Session Info**: Document includes R version and package versions
3. **Exact Methods**: Detailed parameter specifications throughout
4. **Raw Data**: QIIME2 artifacts preserved for reanalysis
5. **Metadata**: Complete sample metadata with collection information

### Citation
If you use this code or data, please cite:

```
Liberman, R., Lopez, J.V., & Coffroth, M.A. (2025). 
Host-specific bacterial microbiomes in four octocoral species from southeastern Florida. 
[Journal name pending]. DOI: [pending]
```

### Session Information
Analysis performed with:
- R version 4.3.0 (or higher)
- microeco version 1.0.0 (or higher)
- Key package versions documented in HTML output

## Project Structure

### Data Flow
```
QIIME2 Processing
    ↓
Raw QIIME2 Artifacts (.qza files)
    ↓
R Analysis Pipeline (this repository)
    ↓
Statistical Outputs + Visualizations
    ↓
Manuscript Figures + Supplementary Data
```

### Analysis Modules

1. **Module 1: Data Import & QC** (Lines 52-230)
   - Import QIIME2 data
   - Filter samples and taxa
   - Generate sequencing statistics

2. **Module 2: Diversity Analysis** (Lines 232-480)
   - Alpha diversity metrics
   - Beta diversity ordination
   - Statistical testing

3. **Module 3: Composition** (Lines 482-640)
   - Taxonomic bar plots
   - Abundance patterns
   - Community structure

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
