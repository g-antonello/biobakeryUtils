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
    
    permanova_slot <- metadata(tse)[[name]]
    
    # get PERMANOVA data frame, this is always generated
    permanova_res.df <- data.frame(
      term = rownames(permanova_slot$permanova),
      PERMANOVA_df = permanova_slot$permanova$Df,
      PERMANOVA_SumOfSqs = permanova_slot$permanova$SumOfSqs,
      PERMANOVA_R2 = permanova_slot$permanova$R2,
      PERMANOVA_P = permanova_slot$permanova[["Pr(>F)"]])
      
    
    if(any(grepl("homogeneity", names(permanova_slot)))){
      # case 1 - homogeneity was tested
        homogeneity_res.df <- 
        data.frame(
          term = rownames(permanova_slot$homogeneity),
          HOMOGENEITY_P = permanova_slot$homogeneity$`Explained variance`,
          HOMOGENEITY_TotVariance = permanova_slot$homogeneity$`Total variance`,
          HOMOGENEITY_ExplVariance = permanova_slot$homogeneity$`Explained variance`
        )
        
      results_clean <- merge(permanova_res.df, homogeneity_res.df, by = "term", all.x = TRUE)
      # keep row order of permanova results
      results_clean <- results_clean[match(permanova_res.df$term, results_clean$term),]
      
      
    } else {
      # case 2, homogeneity was not tested
      results_clean <- permanova_res.df
    }
    
  return(results_clean)
  
}
