## Post-marketing Assessment of Antibody-Drug Conjugates: Proof-of-concept using Trastuzumab-Drug Conjugates, Model-Based Meta-Analysis, and a Clinical Utility Index Approach.

### Overview
This repository contains the study materials and analysis for the project titled "Post-marketing Assessment of Antibody-Drug Conjugates: Proof-of-concept using Trastuzumab-Drug Conjugates, Model-Based Meta-Analysis, and a Clinical Utility Index Approach." This study evaluates the benefit-risk profiles of two approved trastuzumab-drug conjugates, trastuzumab emtansine (T-DM1) and trastuzumab deruxtecan (T-DXd).

### Contents
1. Data Files: Data used in the analysis ("Supplementary_Tables.xlsx" and "base.txt").
2. Scripts: R scripts for data processing, model-based meta-analysis (MBMA), and clinical utility index (CUI) calculations.

### Analysis
Before analysis, install the following R packages:
   ```{r install required packages, include = FALSE}
   # List of required packages
    required_packages <- c("bslib", "ComplexUpset", "export", "ggbump", "mgcv", "patchwork", "PKNCA", "RColorBrewer", "readxl", "shiny", "shinyjs", "shinyWidgets", "tidyverse")

  # Function to install missing packages
  install_if_missing <- function(packages) {
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
      }
    }
  }

  # Install missing packages
  install_if_missing(required_packages)
  ```

The lixoftConnectors also need installation as described at https://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/.  
`r install.packages("C:/Program Files/Lixoft/MonolixSuite2024R1/connectors/lixoftConnectors.tar.gz", repos = NULL, type="source", INSTALL_opts ="--no-multiarch")`

For downstream analysis, six R scripts are provided, but only four need to be run.
* '1. PopPK analysis.R'
* '2. Simulating with uncertainty.R'
* '4. CUI analysis.R' - this script will also source and run '3. ER analysis.R'
* '5b. CUI_shiny.R' - this script will also source and run '5a. Process Shiny data.R'

### Key Findings
Demonstrates the use of MBMA and CUI to assess the benefit-risk profiles of ADCs.
Shows how systematic review and composite scoring can optimize dosing and maximize patient benefit.
Provides insights into the efficacy and safety of T-DM1 and T-DXd for different cancer types.

### Usage
This repository is intended for researchers and practitioners in oncology, pharmacology, and data science. The provided scripts and data can be used to replicate the study's findings or to apply the methodology to other datasets.

### Citation
If you use this repository in your research, please cite the original study as follows:
[Author Names], "Post-marketing Assessment of Antibody-Drug Conjugates: Proof-of-concept using Trastuzumab-Drug Conjugates, Model-Based Meta-Analysis, and a Clinical Utility Index Approach," [Journal Name], [Year].
