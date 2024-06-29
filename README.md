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
Here, we provide three examples for using CompassR from single sample level, cell type level, and multi-sample level.

1. **single-sample level**
Please go the the [vignettes](./vignettes/sample_example.Rmd) for single-sample level example.
2. **multi-sample level**
Please go the the [vignettes](./vignettes/tissue_example.Rmd) for multi-sample level example.
1. **celltype level**
Please go the the [vignettes](./vignettes/celltype_example.Rmd) for celltype level example.
