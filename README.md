# CompassR

## Overview
We developed CompassR to support more in-depth analyses and visualizations of CRE-gene linkages. CompassR is an open source R package built upon [CompassDB](http://compass-db.com/). CompassR can be used to analyze and visualize CRE-gene linkages in one or multiple samples, and identify TFs whose motifs are enriched in selected CREs.

## Installation

To install the latest release of CompassR from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("changxinw/CompassR")
```

## Tutorial
Here, we provide three examples for using CompassR from tissue level, celltype level, and single-sample level.

1. **Tissue level:**
Please go the the [Tissue example](https://changxinw.github.io/CompassR/articles/tissue_example.html) for tissue-level example.
2. **Celltype level:**
Please go the the [Celltype example](https://changxinw.github.io/CompassR/articles/celltype_example.html) for celltype-level example.
3. **Sample level:**
Please go the the [Sample example](https://changxinw.github.io/CompassR/articles/sample_example.html) for sample-level example.
4. **Integrate your inhosue data:**
Learn how to integrat your inhouse data from [Integrate inhouse data](https://changxinw.github.io/CompassR/articles/integrat_inhouse_data.html)
