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

## -----------------------------------------------------------------------------
gene = "HEY2"
celltype = "Cardiomyocyte"
expr_vec = query_exprssion("hg38", gene, celltype = TRUE)
link_res = query_linkage("hg38", gene, celltype = TRUE)

link_df = link_res[["linkage"]]
sample_df = link_res[["samples"]]

sample_df = sample_df[sample_df$bio_source %in% c("Heart right ventricle", "Left cardiac atrium"), ]
sample_df = sample_df[, c("sample_id", "bio_source")]
colnames(sample_df) = c("Sample", "Group")

link_df = link_df[link_df$sample %in% sample_df$Sample & link_df$celltype == celltype, ]
# link_df[, c("sample", "gene", "peak", "zscore", "source", "celltype")]

sample_df = sample_df[sample_df$Sample %in% link_df$sample, ]
sample_df = sample_df[order(sample_df$Group), ]

expr_df = data.frame(t(expr_vec)[paste0(sample_df$Sample, ".", celltype), ])
colnames(expr_df) = c("Gene")
rownames(expr_df) = sapply(rownames(expr_df), function(x) {unlist(strsplit(x, "\\."))[1]})
expr_df$Sample = rownames(expr_df)

link_df$sample_id = link_df$sample
nrow(link_df)

## -----------------------------------------------------------------------------
scheme <- getScheme()
scheme$GdObject$cex.axis = 1
addScheme(scheme, "myScheme")
options(Gviz.scheme="myScheme")

sample_df$Group = c(rep("Ventricle", 2), rep("Atrium", 2))
pdf("Heatmap_celltype_cadiomyocytes_HEY2.pdf", width = 10, height = 6)
genome_track_map(link_df, sample_df, gene, expr_df, assembly = "hg38", legend.position="left", t = -30, b = 30)
dev.off()

