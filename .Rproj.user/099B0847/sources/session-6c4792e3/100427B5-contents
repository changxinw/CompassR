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

#' functions to plot out the genome track
#'
#' @param df gene-cre linkage dataframe
#' @param gene gene to plot
#' @param assembly either mm10 or hg38
#' @param sample_ids samples that want to plot
#' @param min_x
#' @param max_x
#'
#' @return
#' @export
#'
#' @examples
plot_genome_track <- function(df, gene, assembly, sample_ids, min_x = NULL, max_x = NULL){
  bm_host = "https://nov2020.archive.ensembl.org"
  if (assembly == "mm10"){
    bm <- useEnsembl(host = bm_host,
                     biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "mmusculus_gene_ensembl") #hsapiens_gene_ensembl
    grtrack <- BiomartGeneRegionTrack(genome = "mm10",
                                      name = "Gene",
                                      symbol = gene,
                                      biomart = bm,
                                      transcriptAnnotation = "symbol",
                                      collapseTranscripts = "longest",
                                      filters=list(mgi_symbol=gene) ### 'hgnc_symbol''mgi_symbol''uniprot_gn_symbol'
                                      # stacking = "hide"
    )
  } else if (assembly == "hg38"){
    bm <- useEnsembl(host = bm_host,
                     biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl") #hsapiens_gene_ensembl
    grtrack <- BiomartGeneRegionTrack(genome = "hg38",
                                      name = "Gene",
                                      symbol = gene,
                                      biomart = bm,
                                      transcriptAnnotation = "symbol",
                                      collapseTranscripts = "longest",
                                      filters=list(hgnc_symbol=gene) ### 'hgnc_symbol''mgi_symbol''uniprot_gn_symbol'
                                      # stacking = "hide"
    )
  }
  df_gr = df %>%
    pivot_wider(names_from = sample_id, values_from = zscore, id_cols = c('Chromosome', 'Start', 'End')) %>%
    makeGRangesFromDataFrame(seqnames.field="Chromosome", start.field="Start", end.field="End", keep.extra.columns=TRUE, ignore.strand=TRUE)
  df_gr@elementMetadata = df_gr@elementMetadata[, sample_ids]
  genome(df_gr) = assembly
  chr <- as.character(unique(seqnames(df_gr)))
  colors = colorRampPalette(brewer.pal(n = 9, name = "Blues")[5:9])(100)
  cre_track = DataTrack(df_gr, name = "Link", type = "heatmap", showSampleNames = TRUE, cex.sampleNames = 0.6, col.sampleNames = "black", gradient = colors)#, groups = rep(c("male", "female"), each = 4)
  if (is.null(min_x) & is.null(max_x)) {
    min_x = min(min(start(cre_track@range)), min(start(grtrack@range)))
    max_x = max(max(end(cre_track@range)), max(end(grtrack@range)))
    xwidth = max_x - min_x
  } else {
    xwidth = 0
  }
  itrack <- IdeogramTrack(genome = assembly, chromosome = chr)
  gtrack <- GenomeAxisTrack()
  ### add additional track as peaks panel, should combine individual peaks first
  peak_gr = GenomicRanges::reduce(df_gr, min.gapwidth=1L)
  peak_gr$name = 1:length(peak_gr)
  peak_track = AnnotationTrack(peak_gr, id = peak_gr$name, name = "Peak", showFeatureId=TRUE, stacking = "squish", fontcolor.item="black", cex.feature = 0.5)
  # stacking(peak_track) <- "full"
  # p_peak = as.ggplot(~Gviz::plotTracks(peak_track, from = min_x - 0.05*xwidth, to = max_x + 0.05*xwidth, stacking = "dense"))
  ### get plot result
  track_list = list(itrack, gtrack, cre_track, peak_track, grtrack)
  all_tracks = Gviz::plotTracks(track_list, from = min_x - 0.05*xwidth, to = max_x + 0.05*xwidth)
  # gviz_obj = ~Gviz::plotTracks(track_list, from = min_x - 0.05*xwidth, to = max_x + 0.05*xwidth)
  # p_all = as.ggplot(~Gviz::plotTracks(track_list, from = min_x - 0.05*xwidth, to = max_x + 0.05*xwidth))
  graph_width = as.data.frame(all_tracks$titles@coords)
  graph_width$width = graph_width$y2 - graph_width$y1
  return(list(track_list, graph_width, min_x - 0.05*xwidth, max_x + 0.05*xwidth, peak_gr))
}


#' Function to plot annotation bar for each sample
#'
#' @param df_group annotation dataframe with two columns: Sample and Group
#' @param sample_order the vector indicate sample order of all samples
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
plot_annotation_bar = function(df_group, sample_order = df_group$Sample, ...){
  df_group$Sample = factor(df_group$Sample, levels = rev(sample_order))
  if (is.numeric(df_group$Group)){
    p = ggplot(df_group, aes(y = Sample, x = 1, fill = Group)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low="blue", high="red") +
      # scale_fill_brewer(palette = 'Set2') +
      theme_nothing()
  } else {
    p = ggplot(df_group, aes(y = Sample, x = 1, fill = Group)) +
      geom_tile(color = "white") +
      scale_fill_brewer(palette = 'Set2') +
      theme_nothing()
  }
  return(p + theme(...))
}


#' Function to plot out expression of gene as barplot
#'
#' @param sample_ids
#' @param gene
#' @param p_group
#' @param sample_folder
#'
#' @return
#' @export
#'
#' @examples
plot_expr_barplot = function(expr_list, group_df, gene, p_group){
  sample_ids = group_df$Sample
  expr_vec = sapply(sample_ids, function(x){
    expr_list[[x]][gene]
  })
  plot_df = data.frame(Gene = expr_vec, Sample = sample_ids)
  plot_df$Sample = factor(plot_df$Sample, levels = rev(plot_df$Sample))
  rownames(plot_df) = plot_df$Sample
  p <- ggplot(data=plot_df, aes(x=Sample, y=Gene)) +
    labs(y = gene) +
    geom_bar(stat="identity") +
    scale_y_continuous(position = "right", expand = c(0,0)) +
    coord_flip() +
    ylim2(p_group) +
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      axis.text.y=element_blank(),
      axis.title.y=element_blank()
    )
  return(p)
}

#' Function to plot out expression of gene as boxplot
#'
#' @param expr_list_all
#' @param sample_ids
#' @param gene
#' @param p_group
#'
#' @return
#' @export
#'
#' @examples
plot_expr_boxplot = function(expr_list_all, group_df, gene, p_group){
  sample_ids = group_df$Sample
  expr_list = lapply(sample_ids, function(x){
    expr_df = data.frame(expr_list_all[[x]][gene, ])
    expr_df$Sample = x
    colnames(expr_df) = c("Gene", "Sample")
    return(expr_df)
  })
  plot_df = do.call("rbind", expr_list)
  plot_df$Sample = factor(plot_df$Sample, levels = rev(sample_ids))
  p <- ggplot(data=plot_df, aes(x=Sample, y=Gene)) +
    labs(y = gene) +
    # geom_bar(stat="identity") +
    geom_boxplot(fill='#A4A4A4') +
    scale_y_continuous(position = "right", expand = c(0,0)) +
    coord_flip() +
    ylim2(p_group) +
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      axis.text.y=element_blank(),
      axis.title.y=element_blank()
    )
  return(p)
}

#' Title
#'
#' @param track_plot
#' @param annot_plot
#' @param expr_plot
#' @param track_width
#' @param legend.position
#' @param t
#' @param b
#'
#' @return
#' @export
#'
#' @examples
combind_plots = function(track_plot, annot_plot, expr_plot, track_width, legend.position="bottom", t = -20, b = 20){
  group_legend = plot_grid(get_legend(annot_plot + theme(legend.position = legend.position)))
  annot_expr = annot_plot +
    expr_plot +
    plot_layout(widths = c(1, 20)) +
    theme(plot.margin = margin(t = track_width$width[2]+t, b = b))
  design <- "
      1#
      12
      12
      13
      13
    "
  p = track_plot + annot_expr + group_legend + plot_layout(design = design, heights = track_width$width, widths = c(50, 20))
  return(p)
}

#' Function to make the genome track plot
#'
#' @param link_df
#' @param group_df
#' @param gene
#' @param expr_list
#' @param output_file
#' @param assembly
#' @param plot_type
#' @param width
#' @param height
#' @param ... additional parameters passed to plots combination
#'
#' @return
#' @export
#'
#' @examples
genome_track_map = function(link_df, group_df, gene, expr_list, output_file = "./test.pdf", assembly = "hg38", plot_type = "barplot", width = 12, height = 6, ...){
  group_df = group_df[order(group_df$Group), ]
  link_gene_df = link_df[link_df$gene == gene & link_df$sample_id %in% group_df$Sample, ]
  ### prepare dataframe
  if (length(unique(link_gene_df$sample_id))<=2){ return(NULL) }
  ### plot
  annot_plot = plot_annotation_bar(group_df)
  expr_plot = plot_expr_barplot(expr_list, group_df, gene, annot_plot)
  track_plot = plot_genome_track(link_gene_df, gene, assembly, group_df$Sample)
  p_track = ggplotify::as.ggplot(~Gviz::plotTracks(track_plot[[1]], from = track_plot[[3]], to = track_plot[[4]]), envir=environment())
  pdf(output_file, width = width, height = height)
  p = combind_plots(p_track, annot_plot, expr_plot, track_plot[[2]], ...)
  print(p)
  dev.off()
  return(p)
}
