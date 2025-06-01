#' Download and incorporate a cleaned mpa tree version into a Tree Summarized Experiment
#'
#' The function supports the "UNCLASSIFIED" reads option by adding a fake 
#' tip to the root with distance = 1
#'
#' @param data.tse \code{TreeSummarizedExperiment} object without a `rowTree`
#' @param CHOCOPhlAn_version \code{character} that specifies with tree should be 
#' downloaded. supported so far: 202403 (default) and 202307. See other timestamps in http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/
#'
#' @importFrom ape read.tree keep.tip
#' @importFrom TreeTools AddTip
#'
#' @returns A \code{TreeSummarizedExperiment} now including a `rowTree`
#' @export
#'
#' @examples
#' 

AddPhyloTree_to_mpa_tse <- function(data.tse, CHOCOPhlAn_version = "202403"){
  if(!is.null(data.tse)){
    message("Tree already present, returning untouched input")
    return(data.tse)
  }
  trees_available <- c("mpa_vJan25_CHOCOPhlAnSGB_202503.nwk", "mpa_vJun23_CHOCOPhlAnSGB_202403.nwk", "mpa_vJun23_CHOCOPhlAnSGB_202307.nwk")
  
  wanted_tree <- trees_available[grepl(CHOCOPhlAn_version, trees_available)]
  
  if(isFALSE(wanted_tree)){
      stop(paste("timestamps supported are: ", paste(timestamps_available, collapse = ", ")))
  }
  
  mpa.tre <- ape::read.tree(paste0("http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/", trees_available[grepl(CHOCOPhlAn_version, trees_available)]))
  
  #' rename tips as SGBs. Referring to Aitor Blanco's note in the following biobakery
  #' topic: https://forum.biobakery.org/t/inquiry-regarding-metaphlan-sgbs-phylogenetic-tree/4442/3
  mpa.tre$tip.label <- paste0("t__SGB", mpa.tre$tip.label)
  
  # add UNCLASSIFIED if necessary
  if("UNCLASSIFIED" %in% rownames(data.tse)){
    mpa.tre <- AddTip(mpa.tre, label = "UNCLASSIFIED", where = 0)
  }
  
  # subset tree only with wanted tips
  relevant_tips <- intersect(mpa.tre$tip.label, rownames(data.tse))
  tree_subset <- ape::keep.tip(mpa.tre, relevant_tips)
  
  #' maybe not necessary:reorder the TreeSummarized experiment based on tip 
  #' labels order
  data_reordered.tse <- data.tse[tree_subset$tip.label,]
  
  # add tree to the TreeSummarizedExperiment object
  
  rowTree(data_reordered.tse) <- tree_subset
  
  return(data_reordered.tse)
}

