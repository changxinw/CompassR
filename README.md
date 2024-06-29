# CompassR

## Overview
We developed CompassR to support more in-depth analyses and visualizations of CRE-gene linkages. CompassR is an open source R package built upon [CompassDB](http://compass-db.com/). CompassR can be used to analyze and visualize CRE-gene linkages in one or multiple samples, and identify TFs whose motifs are enriched in selected CREs.

## Installation

To install the latest release of Signac from CRAN:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("changxinw/CompassR")
```

## Tutorial
Here, we provide three examples for using CompassR from tissue level, celltype level, and single-sample level.

1. **Tissue level:**
Please go the the [tissue_example.Rmd](./vignettes/tissue_example.Rmd) for tissue-level example.
2. **Celltype level:**
Please go the the [celltype_example.Rmd](./vignettes/celltype_example.Rmd) for celltype-level example.
3. **Sample level:**
Please go the the [sample_example.Rmd](./vignettes/sample_example.Rmd) for sample-level example.
