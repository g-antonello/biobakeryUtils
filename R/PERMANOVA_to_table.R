#' PERMANOVA to table
#'
#' A utility function that produces a data.frame with PERMANOVA and betadisper
#' statistics from a TreeSummarizedExperiment object
#'
#' @param tse \code{TreeSummarizedExperiment} object that has been already populated with `mia::addPERMANOVA`
#' @param name \code{character} specifying the 'name' paramenter that was passed to `mia::addPERMANOVA`
#'
#' @importFrom broom tidy
#' @importFrom dplyr rename left_join transmute
#' @importFrom S4Vectors metadata
#' 
#' @returns A \code{data.frame} with permanova and betadisper + permutest statistics
#' ready to be reported or saved as tables.
#'  
#' @export
#'
#' @examples
#' library(mia)
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Apply relative transformation
#' tse <- transformAssay(tse, method = "relabundance")
#' # Perform PERMANOVA
#' tse <- addPERMANOVA(
#'   tse,
#'   assay.type = "relabundance",
#'   method = "bray",  
#'   formula = x ~ SampleType,
#'   permutations = 99
#' )
#' # The results are stored to metadata
#' metadata(tse)[["permanova"]]
#' 
#' # PERMANOVA_to_table makes it faster to retrieve all pieces of it in the
#' # tse object
#' permanova.df <- PERMANOVA_to_table(tse, name = "permanova")
#' permanova.df


PERMANOVA_to_table <- function(tse, name = "permanova") {
  if (!all(names(metadata(tse)[[name]]) %in% c("permanova", "homogeneity"))) {
    results_table.df <- tidy(metadata(tse)[[name]])
  } else {
    permanova.df <- suppressWarnings(tidy(metadata(tse)[[name]][["permanova"]]))
    permanova.df <- rename(permanova.df, p.value_PERMANOVA = p.value)
    
    tmp <- as.data.frame(metadata(tse)[[name]][["homogeneity"]])
    betadisper.df <- transmute(tmp,
        term = rownames(tmp),
        Tot.variance_betadisper = `Total variance`,
        Expl.variance_betadisper = `Explained variance`,
        p.value_betadisper = `Pr(>F)`
      )
    
    results_table.df <- left_join(permanova.df, betadisper.df, by = "term")
  }
  
  return(results_table.df)
  
}