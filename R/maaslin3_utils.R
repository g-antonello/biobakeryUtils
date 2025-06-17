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
#' @importFrom dplyr full_join arrange
#' @importFrom SummarizedExperiment rowData
#' 
#' @export
#'
#' @examples
#' 
#' library(maaslin3)
#' WallenZD_2022.tse <- mia::importMetaPhlAn(
#' file = system.file("extdata",
#'                    "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                    package = "biobakeryUtils"),
#' col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", package = "biobakeryUtils"))
#'
#' maaslin3_results_raw <- maaslin3(
#'   input_data = WallenZD_2022.tse, 
#'   formula = ~ Case_status + Sex + Age_at_collection, 
#'   transform = "LOG", 
#'   output = tempdir(),
#'   normalization = "TSS", # scale values per samples from 0 to 1
#'   min_prevalence = 0.05, 
#'   min_abundance = 0.001, 
#'   standardize = FALSE, 
#'   verbosity = "ERROR" 
#' )
#' 
#' maaslin3_with_taxonomy <- add_taxonomy_to_maaslin3(maaslin3_res = maaslin3_results_raw, input.tse = WallenZD_2022.tse, taxLevel = "Guess")
#' 
#' head(maaslin3_with_taxonomy$abundance)


add_taxonomy_to_maaslin3 <- function(maaslin3_res, input.tse, taxLevel = "Guess"){
  
  # prepare maaslin3 data
    results.list <- list(
      abundance = maaslin3_res$fit_data_abundance$results,
      prevalence = maaslin3_res$fit_data_prevalence$results
    )

  # if taxLevel is not explicitly stated, guess it
  if(taxLevel == "Guess"){
    
    prefix <- unique(grep(paste(names(all_taxonomy_levels), collapse = "|"),substr(results.list$abundance$feature, 1, 3), value = TRUE))
   
    if(length(prefix) != 1){
      stop("Taxonomy is not guessable")
    }
    taxLevel <- all_taxonomy_levels[prefix]
  }
  
  # get taxonomic information up to the taxonomic level needed
    taxonomy_info <- as.data.frame(rowData(input.tse))
    taxonomy_info_subset <- taxonomy_info[,1:which(agrepl(taxLevel, colnames(taxonomy_info)))]
    colnames(taxonomy_info_subset)[ncol(taxonomy_info_subset)] <- "feature"
    
  # merge results with taxonomy, keeping them split
  results_w_taxonomy.list <- lapply(results.list, function(x) {
      
      return(
        arrange(
          full_join(
            taxonomy_info_subset, x, by = "feature"),
          qval_individual)
        )
    }
    )
  
  return(results_w_taxonomy.list)
}



#' Write curated maaslin3 output
#'
#' This is a quick function to save data that result from the previous function
#' called `add_taxonomy_to_maaslin3`
#' 
#' @param maaslin3_res_list \code{list} of results as come out of 
#' `add_taxonomy_to_maaslin3`
#' @param out.dir \code{character} specifying where to save results
#'  
#' @importFrom readr write_tsv
#' @returns NULL, nothing is returned. Files are save as 
#' `all_results_w_taxonomy_abundance.tsv` and 
#' `all_results_w_taxonomy_prevalence.tsv`
#' 
#' @export
#'
#' @examples
#'
#' library(maaslin3)
#' WallenZD_2022.tse <- mia::importMetaPhlAn(
#' file = system.file("extdata",
#'                    "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                    package = "biobakeryUtils"),
#' col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", package = "biobakeryUtils"))
#'
#' maaslin3_results_raw <- maaslin3(
#'   input_data = WallenZD_2022.tse, 
#'   formula = ~ Case_status + Sex + Age_at_collection, 
#'   transform = "LOG", 
#'   output = tempdir(),
#'   normalization = "TSS", # scale values per samples from 0 to 1
#'   min_prevalence = 0.05, 
#'   min_abundance = 0.001, 
#'   standardize = FALSE, 
#'   verbosity = "ERROR" 
#' )
#' 
#' maaslin3_with_taxonomy <- add_taxonomy_to_maaslin3(maaslin3_res = maaslin3_results_raw, input.tse = WallenZD_2022.tse, taxLevel = "Guess")
#' 
#' head(maaslin3_with_taxonomy$abundance)
#' head(maaslin3_with_taxonomy$prevalence)
#' 
#' write_maaslin3_curated_tables(maaslin3_with_taxonomy, tempdir())
#' 
#' list.files(tempdir(), pattern = "all_results_w_taxonomy*")

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
