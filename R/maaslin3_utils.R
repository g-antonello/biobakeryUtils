#' Add taxonomy to maaslin3 results
#' 
#' This function adds taxonomies to raw maaslin3 results by using the original
#' input `TreeSummarizedExperiment` microbiome input data. So far, it is designed
#' to work only with MetaPhlAn4 output tables
#'
#' @param maaslin3_res \code{maaslin3 raw output list} after running the model.
#' 
#' @param input.tse The \code{TreeSummarizedExperiment} microbiome input data 
#' used for modeling
#' 
#' @param taxLevel \code{character} specifying the taxonomic level. Default is 
#' "Guess", which tries to guess based on the taxa prefixes (eg: Species = s__).
#'
#' @returns A \code{list} of \code{data.frame}s, one for the abundance and one 
#' for the prevalence results, with taxonomy in the first columns and maaslin3 
#' results following. Taxa not tested are also included for completeness.
#' 
#' @export
#'
#' @examples
#' 

add_taxonomy_to_maaslin3 <- function(maaslin3_res, input.tse, taxLevel = "Guess"){
  
  # prepare maaslin3 data
    results.list <- list(
      abundance = maaslin3_res$fit_data_abundance$results,
      prevalence = maaslin3_res$fit_data_prevalence$results
    )

  # if taxLevel is not explicitly stated, guess it
  if(taxLevel == "Guess"){
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
    
    prefix <- unique(grep(paste(names(all_possible_taxonomy_colnames), collapse = "|"),substr(results.list$abundance$feature, 1, 3), value = TRUE))
   
    if(length(prefix) != 1){
      stop("Taxonomy is not guessable")
    }
    taxLevel <- all_possible_taxonomy_colnames[prefix]
  }
  
  # get taxonomic information up to the taxonomic level needed
    taxonomy_info <- as.data.frame(rowData(input.tse))
    taxonomy_info_subset <- taxonomy_info[,1:which(colnames(taxonomy_info) == taxLevel)]
    colnames(taxonomy_info_subset)[ncol(taxonomy_info_subset)] <- "feature"
    
  # merge results with taxonomy, keeping them split
  results_w_taxonomy.list <- lapply(results.list, function(x) {
      
      return(
        dplyr::full_join(taxonomy_info_subset, x, by = "feature") %>%
          arrange(qval_individual))}
    )
  
  return(results_w_taxonomy.list)
}



#' Write curated maaslin3 output
#'
#' @param maaslin3_res_list \code{list} of results as come out of 
#' `add_taxonomy_to_maaslin3`
#' @param out.dir \code{character} specifying where to save results
#'
#' @returns NULL, nothing is returned. Files are save as .tsv
#' @export
#'
#' @examples
#' 
#' 

write_maaslin3_curated_tables <- function(maaslin3_res_list, out.dir){
  
  write_tsv(
    maaslin3_res_list$abundance,
    file = file.path(out.dir, "all_results_w_taxonomy_abundance.tsv")
  )
  
  write_tsv(
    maaslin3_res_list$abundance,
    file = file.path(out.dir, "all_results_w_taxonomy_prevalence.tsv")
  )
}
