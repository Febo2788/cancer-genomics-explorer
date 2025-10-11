# Cancer Genomics Explorer

An interactive Shiny web application for exploring cancer genomics data from The Cancer Genome Atlas (TCGA). Visualize gene expression patterns, perform survival analysis, and discover gene correlations across multiple cancer types.

## [Live Demo](https://febo2788.shinyapps.io/cancer_genomics_app/)

![R](https://img.shields.io/badge/R-4.0%2B-blue)
![Shiny](https://img.shields.io/badge/Shiny-1.7%2B-brightgreen)

## Features

### 1. Gene Expression Visualization
- Compare gene expression levels across different cancer types
- Toggle between boxplot and violin plot visualizations
- Interactive sample-level data points
- Summary statistics for selected genes

### 2. Survival Analysis
- Kaplan-Meier survival curves comparing high vs. low gene expression
- Adjustable expression cutoff (percentile-based)
- Log-rank test p-values
- Risk tables and confidence intervals
- Median survival time comparisons

### 3. Gene Correlation Analysis
- Identify genes correlated with your gene of interest
- Top 20 most correlated genes visualization
- Pearson correlation coefficients
- Distinguish positive and negative correlations

### 4. Interactive Data Exploration
- Three cancer types: Breast (BRCA), Lung Adenocarcinoma (LUAD), Colon (COAD)
- 25 common cancer genes (TP53, KRAS, EGFR, BRCA1, etc.)
- Downloadable high-resolution plots
- Real-time responsive visualizations

## Installation

### Prerequisites
- R version 4.0 or higher
- RStudio (recommended)

### Required R Packages

```r
# CRAN packages
install.packages(c("shiny", "ggplot2", "dplyr", "survival", "survminer", "tidyr"))

# Bioconductor packages (for TCGA data download)
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
```

## Quick Start

1. **Clone or download this repository**

```bash
git clone <repository-url>
cd cancer_genomics_app
```

2. **Download TCGA data**

```r
setwd("path/to/cancer_genomics_app")
source("download_quick_fixed.R")
```

This downloads real TCGA data (50 samples per cancer type) and creates `demo_data.rds`.

3. **Launch the app**

```r
library(shiny)
runApp()
```

The app will open in your default web browser.

## Usage Guide

### Selecting Data
1. Choose a **Cancer Type** from the dropdown (BRCA, LUAD, or COAD)
2. Select a **Gene** to analyze from the list of common cancer genes

### Expression Tab
- View gene expression distributions across all cancer types
- Switch between boxplot and violin plot styles
- Read summary statistics in the text output

### Survival Tab
- Analyze patient survival based on gene expression levels
- Adjust the cutoff percentile to redefine high/low expression groups
- Interpret the Kaplan-Meier curves and p-value

### Correlations Tab
- Discover which genes are most correlated with your selected gene
- Positive correlations (green) vs. negative correlations (red)
- View detailed correlation coefficients in the table

### Downloading Plots
- Each tab has a "Download Plot" button
- Saves high-resolution PNG images (300 DPI)
- Suitable for presentations and publications

## Data Source

This application uses data from **The Cancer Genome Atlas (TCGA)**, a landmark cancer genomics program that molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types.

- **Expression values**: log2(FPKM + 1) normalized gene expression
- **Clinical data**: Survival times, vital status, age, gender, stage
- **Cancer types included**:
  - BRCA: Breast invasive carcinoma
  - LUAD: Lung adenocarcinoma
  - COAD: Colon adenocarcinoma

Learn more: [https://www.cancer.gov/tcga](https://www.cancer.gov/tcga)

## Technical Details

### Architecture
- **UI**: Fluid page layout with sidebar controls and tabbed main panel
- **Server**: Reactive programming for real-time data updates
- **Plotting**: ggplot2 for expression/correlation, survminer for survival curves

### Key Functions
- `survfit()`: Kaplan-Meier survival curve fitting
- `Surv()`: Survival object creation
- `cor()`: Pearson correlation calculation
- `ggsurvplot()`: Enhanced survival plot visualization

### File Structure
```
cancer_genomics_app/
├── app.R                      # Main Shiny application
├── download_quick_fixed.R     # TCGA data downloader
├── demo_data.rds              # TCGA dataset (generated)
├── README.md                  # This file
└── .gitignore                 # Git ignore rules
```

## Troubleshooting

**Error: Data file not found**
- Run `source("download_quick_fixed.R")` to download TCGA data

**Plots not showing**
- Ensure all required packages are installed
- Check that demo_data.rds exists in the app directory

**Survival plot shows "Insufficient data"**
- This occurs when clinical data is incomplete for certain samples
- Try a different cancer type or gene

## Future Enhancements

Potential features for future versions:
- Additional cancer types (PRAD, LIHC, SKCM, etc.)
- Genome-wide gene selection (all ~20,000 genes)
- Differential expression analysis
- Pathway enrichment analysis
- Mutation data integration

## License

This project is open source and available for educational purposes.

## Citation

If you use this app in your research, please cite TCGA:

> The Cancer Genome Atlas Research Network. (2013). The Cancer Genome Atlas Pan-Cancer analysis project. *Nature Genetics*, 45(10), 1113-1120.

---

**Built with R Shiny** | **Data from TCGA** | **For Research & Education**
