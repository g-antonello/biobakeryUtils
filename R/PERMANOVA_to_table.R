#' PERMANOVA to table
#'
#' A utility function that produces a data.frame with PERMANOVA and betadisper
#' statistics from a TreeSummarizedExperiment object
#'
#' @param tse \code{TreeSummarizedExperiment} object that has been already 
#' populated with `mia::addPERMANOVA`
#' @param name \code{character} specifying the `name` paramenter that was passed 
#'to `mia::addPERMANOVA` as name to store the test in the `metadata` slot
#' 
#' @returns A \code{data.frame} with permanova and, if tested in addPERMANOVA, 
#' betadisper + permutest statistics ready to be reported or saved as tables. 
#'
#' @importFrom S4Vectors metadata
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
  
    results_table.df <- as.data.frame(metadata(tse)[[name]])
    
    # case 1 - homogeneity was tested
    if(any(grepl("homogeneity", colnames(results_table.df)))){
      results_clean <- data.frame(
        term = rownames(results_table.df),
        df = results_table.df$permanova.Df,
        SumOfSqs = results_table.df$permanova.SumOfSqs,
        R2_stat = results_table.df$permanova.R2,
        p.value_PERMANOVA = results_table.df$permanova.Pr..F.,
        Tot.variance_betadisper = results_table.df$homogeneity.Total.variance,
        Expl.variance_betadisper = results_table.df$homogeneity.Explained.variance,
        p.value_betadisper = results_table.df$homogeneity.Pr..F.
      )  
      
      # case 2, homogeneity was not tested
    } else {
      results_clean <- data.frame(
        term = rownames(results_table.df),
        df = results_table.df$permanova.Df,
        SumOfSqs = results_table.df$permanova.SumOfSqs,
        R2_stat = results_table.df$permanova.R2,
        p.value_PERMANOVA = results_table.df$permanova.Pr..F.
      )  
    }
    
  return(results_clean)
  
}