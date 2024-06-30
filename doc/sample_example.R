## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CompassR)
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

## -----------------------------------------------------------------------------
### Sample level analysis
output = "./"
options(timeout=600)
obj = readRDS(url("http://compass-db.com/static/seurat_object/all/GSM5065524_GSM5065525.rds", "rb"))
download.file("http://compass-db.com/static/fragments/all/GSM5065524_GSM5065525_atac_fragments.tsv.gz",
              paste0(output, "GSM5065524_GSM5065525_atac_fragments.tsv.gz"))
download.file("http://compass-db.com/static/fragments/all/GSM5065524_GSM5065525_atac_fragments.tsv.gz.tbi",
              paste0(output, "GSM5065524_GSM5065525_atac_fragments.tsv.gz.tbi"))


## -----------------------------------------------------------------------------
DefaultAssay(obj) = "ATAC"
frag = Fragments(obj)[[1]]
frag = UpdatePath(frag, new.path = paste0(output, "GSM5065524_GSM5065525_atac_fragments.tsv.gz"), verbose = TRUE)
obj@assays$ATAC@fragments[[1]] = frag

gene = "CD79B"
p = CoveragePlot(
  object = obj,
  group.by = "annot",
  region = gene,
  features = gene,
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 10000,
  extend.downstream = 10000
)
# ggsave("CoveragePlot_CD79B_sample_eg.pdf", p, width = 6, height = 8)

