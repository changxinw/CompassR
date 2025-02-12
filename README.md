# CompassR

## Overview
We developed CompassR to support more in-depth analyses and visualizations of CRE-gene linkages. CompassR is an open source R package built upon [CompassDB](http://compass-db.com/). CompassR can be used to analyze and visualize CRE-gene linkages in one or multiple samples, and identify TFs whose motifs are enriched in selected CREs.

## Installation

Installation of dependent:

``` r
# List of additional packages to install
cran_pkgs <- c("aplot", "cowplot", "dplyr", "tidyr", "ggplot2", "ggplotify", "patchwork", "RColorBrewer", "jsonlite", "httr")

# Install missing packages
sapply(cran_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
})

# Install Bioconductor packages if needed
bioc_pkgs <- c("biomaRt", "GenomicRanges", "Gviz")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioc_pkgs, update = TRUE, ask = FALSE)
```

To install the latest release of CompassR from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("changxinw/CompassR")
```

## Tutorial
Here, we provide three examples for using CompassR from tissue level, celltype level, and single-sample level.

1. **Tissue level:**
Please go the the [Myh6 Regulation in Heart](https://changxinw.github.io/CompassR/articles/tissue_example.html) for tissue-level example.
2. **Celltype level:**
Please go the the [HEY2 Regulation in Cardiomyocytes](https://changxinw.github.io/CompassR/articles/celltype_example.html) for celltype-level example.
3. **Sample level:**
Please go the the [CD79B Regulation in B cells](https://changxinw.github.io/CompassR/articles/sample_example.html) for sample-level example.
4. **Integrate your inhosue data:**
Learn how to integrat your inhouse data from [Integrating In-house Data with Datasets in CompassDB](https://changxinw.github.io/CompassR/articles/integrat_inhouse_data.html)
