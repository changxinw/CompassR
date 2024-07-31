
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
    summarise(n = n()) %>%
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

#' wrap up function to generate results of examples quickly
#'
#' @param gene gene involved in each example
#' @param output output folder for download files
#'
#' @return figures for example results
#'
#' @examples NA
quick_example_compassr = function(gene, output = "./"){
  if (gene == 'Myh6'){
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
    expr_df = data.frame(t(expr_vec)[sample_df$Sample, ])
    colnames(expr_df) = c("Gene")
    expr_df$Sample = rownames(expr_df)
    link_df$sample_id = link_df$sample
    genome_track_map(link_df, sample_df, "Myh6", expr_df, assembly = "mm10", legend.position="left", t = -20, b = 20)
    track_plot = plot_genome_track(link_df, "Myh6", "mm10", sample_df$Sample)
    peaks = track_plot[[length(track_plot)]]
    peaks$name = as.character(peaks$name)
    peak_order = make_peak_group(link_df, peaks, sample_df)
    ht_peak_order = peak_order[peak_order[, 2] >= 0.5 & peak_order[, 1]<0.5 & peak_order[, 3]<0.5 & peak_order[, 4]<0.5, ]
    ht_peaks = rownames(ht_peak_order)
    peak_vec = sub(":", "-", as.character(peaks))
    ht_peaks = peak_vec[as.numeric(ht_peaks)]
    tf = tf_binding("mm10", paste0(ht_peaks, collapse = "_"))
    plot_giggle(tf)
    # print(p)
  } else if (gene == "HEY2"){
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
    sample_df = sample_df[sample_df$Sample %in% link_df$sample, ]
    sample_df = sample_df[order(sample_df$Group), ]
    expr_df = data.frame(t(expr_vec)[paste0(sample_df$Sample, ".", celltype), ])
    colnames(expr_df) = c("Gene")
    rownames(expr_df) = sapply(rownames(expr_df), function(x) {unlist(strsplit(x, "\\."))[1]})
    expr_df$Sample = rownames(expr_df)
    link_df$sample_id = link_df$sample
    scheme <- getScheme()
    scheme$GdObject$cex.axis = 1
    addScheme(scheme, "myScheme")
    options(Gviz.scheme="myScheme")
    sample_df$Group = c(rep("Ventricle", 2), rep("Atrium", 2))
    genome_track_map(link_df, sample_df, gene, expr_df, assembly = "hg38", legend.position="left", t = -30, b = 30)
    track_plot = plot_genome_track(link_df, "HEY2", "hg38", sample_df$Sample)
    peaks = track_plot[[length(track_plot)]]
    peaks$name = as.character(peaks$name)
    peak_order = make_peak_group(link_df, peaks, sample_df)
    ht_peak_order = peak_order[peak_order[, 1] >= 0.5 & peak_order[, 2]<0.5, ]
    ht_peaks = rownames(ht_peak_order)
    peak_vec = sub(":", "-", as.character(peaks))
    ht_peaks = peak_vec[as.numeric(ht_peaks)]
    tf = tf_binding("hg38", paste0(ht_peaks, collapse = "_"))
    plot_giggle(tf)
  } else if (gene == "CD79B"){
    ### Sample level analysis
    options(timeout=600)
    obj = readRDS(url("http://compass-db.com/static/seurat_object/all/GSM5065524_GSM5065525.rds", "rb"))
    download.file("http://compass-db.com/static/fragments/all/GSM5065524_GSM5065525_atac_fragments.tsv.gz",
                  paste0(output, "GSM5065524_GSM5065525_atac_fragments.tsv.gz"))
    download.file("http://compass-db.com/static/fragments/all/GSM5065524_GSM5065525_atac_fragments.tsv.gz.tbi",
                  paste0(output, "GSM5065524_GSM5065525_atac_fragments.tsv.gz.tbi"))
    DefaultAssay(obj) = "ATAC"
    frag = Fragments(obj)[[1]]
    frag = UpdatePath(frag, new.path = paste0(output, "GSM5065524_GSM5065525_atac_fragments.tsv.gz"), verbose = TRUE)
    obj@assays$ATAC@fragments[[1]] = frag
    gene = "CD79B"
    CoveragePlot(
      object = obj,
      group.by = "annot",
      region = gene,
      features = gene,
      annotation = TRUE,
      peaks = TRUE,
      extend.upstream = 10000,
      extend.downstream = 10000
    )
    link_df_obj = Links(obj)
    gene_link = link_df_obj[link_df_obj$gene == "CD79B"]
    peak_vec = sub(":", "-", gene_link$peak)
    tf = tf_binding("hg38", paste0(peak_vec, collapse = "_"))
    plot_giggle(tf)
  } else {
    return(NULL)
  }
}
