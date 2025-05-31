#' Taxonomy from vector to data.frame
#' 
#' Expand a character vector of taxonomies into a full taxonomy table based on 
#' the separator provided and the prefixes at each level. Tested on metaphlan4
#' output tables
#'
#' @param taxonomy.chr A \code{character} vector with the taxonomy
#' @param sep A \code{character} specifying the separator to use to split `taxonomy.chr`
#' @param row.names Either `'short'`, `'long'`, or `'none'`. Default is 
#' `'short'`, which returns the highest taxonomic resolution the input 
#' `taxonomy.chr` has for that row. `'long'` returns `taxonomy.chr` as column 
#' names. `'none'` returns a data.frame without rown ames
#
#' @returns A \code{data.frame}
#' @export
#'
#' @examples

expand_taxonomy <- function(taxonomy.chr, sep = "|", row.names = "short"){
  taxonomy.list <- strsplit(taxonomy.chr, sep, fixed = TRUE)
  taxonomy.df <- as.data.frame(Reduce(rbind, taxonomy.list))
  taxonomies_per_column <- apply(taxonomy.df, 2, function(x) unique(substr(x[!grepl("^UN*", x)], 1, 3)))
  all_possible_taxonomy_colnames <- c(
    d__ = "Domain",
    k__ = "Kingdom",
    p__ = "Phylum",
    c__ = "Class",
    o__ = "Order",
    f__ = "Family",
    g__ = "Genus",
    s__ = "Species",
    t__ = "SGB"
  )
  final_colnames <- all_possible_taxonomy_colnames[taxonomies_per_column]
  
  colnames(taxonomy.df) <- final_colnames
  
  if(row.names == "short"){
    rownames(taxonomy.df) <- sapply(taxonomy.list, function(x) x[length(x)])
  } 
  
  if(row.names == "long"){
    rownames(taxonomy.df) <- taxonomy.chr
  }
  
  if(row.names == "none"){
    rownames(taxonomy.df) <- NULL
  }
  
  return(taxonomy.df)
}

#' Taxonomy from data.frame to vector
#'
#' @param taxonomy.df \code{data.frame} with taxonomic levels on columns and individual taxa on rows. Taxa levels must be ordered from Kingdom/Domain to make sense, but the function does not check that
#' @param sep \code{character} specifying the separator avoid `c("_", "-", " ", ",", "\t")` separators, because they will likely be present in the taxonomy table already or because it would conflict with saving the taxonomy table in .csv or .tsv format
#'
#' @returns A \code{character} of complete taxonomies separated by `sep`
#' @export
#'
#' @examples

collapse_taxonomy <- function(taxonomy.df, sep = "|"){
  
  invalidSeps <- c("_", "-", " ", ",", "\t")
  
  if(sep %in% invalidSeps) {
    warning(paste0("Please avoid this separator: ", sep))
  }
  
  return(apply(as.data.frame(rowData(project4_reordered.tse)), 1, function(x) paste(x, collapse = "|")))
  
}

#' Enhanced row renaming in a TreeSummarizedExperiment 
#' 
#' This includes renaming not only rows of a TSE, but renames tree tips as well
#'
#' @param data.tse \code{(Tree)SummarizedExperiment} with or without a phylogenetic tree
#' @param new.rownames = A \code{character} vector of new rownames. IMPORTANT! Make sure the order of old and new row names makes sense before running this function
#'
#' @returns
#' @export
#'
#' @examples

rename_rownames.tse <- function(data.tse, new.rownames){
  
  if(!identical(rownames(data.tse), rowTree(data.tse)$tip.label)){
    stop("Your .tse object does not have identical rownames and tip names in identical order. Consider running `reorder_taxa_with_phyloTree_labels` first")
  }
  
  if((class(data.tse) != "TreeSummarizedExperiment") | is.null(rowTree(data.tse))){
    rownames(data.tse) <- new.rownames
    return(data.tse)
  }

  rownames(data.tse) <- new.rownames
  rowTree(data.tse)$tip.label <- new.rownames
  
  return(data.tse)
}

#' Reorder rows and tip labels of a TreeSummarizedExperiment
#'
#' @param data.tse \code{TreeSummarizedExperiment} object. It must contain a phylogenetic tree
#'
#' @returns \code{TreeSummarizedExperiment} object with identical rownames and tip labels order
#' @export
#'
#' @examples

reorder_taxa_with_phyloTree_labels <- function(data.tse) {
  
  if(is.null(rowTree(data.tse))){
    stop("No phyloTree to reorder taxa with")
  }
  
  data_reordered.tse <- data.tse[rowTree(data.tse)$tip.label,]
  return(data_reordered.tse)
}
