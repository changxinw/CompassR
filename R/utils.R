
#' check if a gene is valid
#'
#' @param assembly hg38 or mm10
#' @param gene gene symbol
#'
#' @return whether the gene is valid
#' @export
#'
#' @examples NA
check_gene = function(assembly, gene){
  return(TRUE)
}

#'  Plot motif or tf binding information
#'
#' @param motif_df data frame of motif or tf binding information
#' @param top_n top n sets to be showed
#' @param title title of the plot
#'
#' @return ggplot2 object of motif plot
#' @export
#' @import ggplot2
#'
#' @examples NA
plot_giggle = function(motif_df, top_n = 15, title = ""){
  motif_df = motif_df[order(-motif_df$combo_score), ]
  motif_df$Factor = factor(motif_df$Factor, levels = rev(motif_df$Factor))
  motif_df = motif_df[1:top_n, ]
  p <- ggplot(data=motif_df, aes(x=Factor, y=combo_score)) +
    labs(x = "TFs", y = "Giggle Score", title = title) +
    geom_bar(stat="identity") +
    coord_flip()
  return(p)
}

#' separate peaks based on their CRE-gene linkages on diferent groups
#'
#' @param link_gene_df CRE-gene link data frame
#' @param peaks coordinate of peaks
#' @param df_group group information of each sample
#'
#' @return data frame of the peak proportion in different groups
#' @export
#'
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr group_by summarise %>% n
#'
#' @examples NA
make_peak_group = function(link_gene_df, peaks, df_group){
  link_gr = makeGRangesFromDataFrame(link_gene_df, seqnames.field="Chromosome", start.field="Start", end.field="End", strand.field="strand", keep.extra.columns=TRUE, ignore.strand=FALSE)
  link_gr$reduced_peak = NA
  olp_df = as.data.frame(findOverlaps(link_gr, peaks))
  for (i in 1:nrow(olp_df)){
    link_gr[olp_df[i, 1]]$reduced_peak = sub(":", "-", peaks[olp_df[i, 2]]$name)
  }
  sample_peak = unique(mcols(link_gr)[, c("sample_id", "reduced_peak")])
  sample_peak$Group = df_group[sample_peak$sample_id, "Group"]
  # sample_peak$count = 1
  peak_group = as.data.frame(sample_peak) %>%
    group_by(Group, reduced_peak) %>%
    summarise(n = dplyr::n()) %>%
    pivot_wider(id_cols = "reduced_peak", names_from = "Group", values_from = n) %>%
    as.data.frame()
  rownames(peak_group) = peak_group[, 1]
  peak_group = peak_group[, -1]
  peak_group[is.na(peak_group)] = 0
  group_freq = summary.factor(df_group$Group)
  for (i in names(group_freq)){
    peak_group[, i] = peak_group[, i] / group_freq[i]
  }
  peak_order = peak_group[order(peak_group[, 1], peak_group[, 2], decreasing = TRUE), ]
  return(peak_order)
}
