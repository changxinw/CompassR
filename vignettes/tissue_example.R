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


## Example for Myh6 in mouse heart sample
expr_vec = query_exprssion("mm10", "Myh6")
link_res = query_linkage("mm10", "Myh6")

link_df = link_res[["linkage"]]
sample_df = link_res[["samples"]]

sample_df = sample_df[sample_df$donor_gender == "Male", ]
sample_df = sample_df[, c("sample_id", "bio_source")]
colnames(sample_df) = c("Sample", "Group")
sample_df = sample_df[order(sample_df$Group), ]

link_df = link_df[link_df$sample %in% sample_df$Sample, ]

## -----------------------------------------------------------------------------
expr_df = data.frame(t(expr_vec)[sample_df$Sample, ])
colnames(expr_df) = c("Gene")
expr_df$Sample = rownames(expr_df)

## -----------------------------------------------------------------------------
link_df$sample_id = link_df$sample
genome_track_map(link_df, sample_df, "Myh6", expr_df, assembly = "mm10", legend.position="left", t = -20, b = 20)

