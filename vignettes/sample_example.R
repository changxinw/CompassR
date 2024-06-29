## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scomicsdb)
library(httr)
library(jsonlite)
library(ggplot2)
library(cowplot)
library(ggforce)
library(dplyr)
library(Seurat)
library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggdendro)
library(cowplot)
library(ggtree) # install with `devtools::install_github("YuLab-SMU/ggtree")` as you need a version newer than what bioconductor serves
library(patchwork)
library(aplot)
library(ggplotify)
library(RColorBrewer)
library(Signac)

