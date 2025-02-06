#' query metadata of samples
#'
#' @param api website api
#'
#' @return data frame of sample information
#' @export
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#'
#' @examples NA
query_meta = function(api = "http://compass-db.com/metaapi/"){
  res = GET(api)
  data = fromJSON(rawToChar(res$content))
  sample_df = as.data.frame(data)
  rownames(sample_df) = sample_df[, 1]
  return(sample_df)
}

#' query linkage of a gene or a genome region
#'
#' @param assembly hg38 or mm10
#' @param gene gene symbol
#' @param coord genome coordinate in the form of chr-start-end
#' @param celltype TRUE or FALSE to indicate if this is cell type specific linkage
#'
#' @return a list with linkage and corresponding sample information
#' @export
#'
#' @examples NA
query_linkage = function(assembly, gene="", coord="", celltype = FALSE){
  api = "http://compass-db.com/linkapi/"
  # check gene name
  if (assembly != "hg38" & assembly != "mm10"){
    print("Assembly must choose between hg38 and mm10!")
  }
  if (gene !="") {
    check_gene(assembly, gene)
  }
  if (celltype){
    parames = paste0("?asm=", assembly, "&gene=", gene, "&interval=", coord, "&source=", "&celltype=True")
  } else {
    parames = paste0("?asm=", assembly, "&gene=", gene, "&interval=", coord, "&source=")
  }
  url = paste0(api, parames)
  url_get = GET(url)
  data = fromJSON(rawToChar(url_get$content))
  link_df = as.data.frame(data$data)
  sample_info = query_meta()
  res_list = list()
  res_list[["linkage"]] = link_df
  sample_df = sample_info[unique(link_df$sample), ]
  sample_df = sample_df[!is.na(sample_df$sample_id), ]
  res_list[["samples"]] = sample_df
  return(res_list)
}


#' Get expression of each sample
#'
#' @param assembly hg38 or mm10
#' @param gene gene symbol
#' @param celltype TRUE or FALSE to indicate if this is cell type specific expression
#'
#' @return data frame of expression
#' @export
#'
#' @examples NA
query_exprssion = function(assembly, gene, celltype = FALSE){
  api = "http://compass-db.com/exprapi/"
  # check gene name
  if (assembly != "hg38" & assembly != "mm10"){
    print("Assembly must choose between hg38 and mm10!")
  }
  if (gene !="") {
    check_gene(assembly, gene)
  }
  if (celltype){
    parames = paste0("?asm=", assembly, "&gene=", gene, "&celltype=True")
  } else {
    parames = paste0("?asm=", assembly, "&gene=", gene)
  }
  url = paste0(api, parames)
  url_get = GET(url)
  data = fromJSON(rawToChar(url_get$content))
  expr_df = as.data.frame(data)
  return(expr_df)
}


#' Obtain TF binding information of genome region(s)
#'
#' @param assembly hg38/mm10
#' @param coords genome coordinate in the form of chr-start-end
#'
#' @return a dataframe of transcription factor binding
#' @export
#'
#' @examples NA
tf_binding = function(assembly, coords){
  api = "http://compass-db.com/tfbsapi/"
  # check gene name
  if (assembly != "hg38" & assembly != "mm10"){
    print("Assembly must choose between hg38 and mm10!")
  }
  # if (gene !="") {
  #   check_gene(assembly, gene)
  # }
  coords = paste0(coords, collapse = "_")
  parames = paste0("?asm=", assembly, "&coords=", coords)
  url = paste0(api, parames)
  url_get = GET(url)
  data = fromJSON(rawToChar(url_get$content))
  tfbs_df = as.data.frame(data)
  return(tfbs_df)
}
