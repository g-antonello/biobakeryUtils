#' get Features stats
#'
#' @param tse \code{TreeSummarizedExperiment} to calculate the stats on
#' @param assay.type \code{character} of the assay to calculate the estimates on (default: `NULL`)
#'
#' @returns \code{data.frame} with one row per feature and some common useful 
#' summary stats in columns
#' 
#' @importFrom tibble tibble rownames_to_column
#' @importFrom TreeSummarizedExperiment assay assayNames
#' 
#' @export
#'
#' @examples
#' library(mia)
#' WallenZD_2022.tse <- mia::importMetaPhlAn(
#' file = system.file("extdata",
#'                    "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                    package = "biobakeryUtils"),
#' col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", 
#' package = "biobakeryUtils"))
#' 
#' Wallen_stats.df <- getFeatureStats(WallenZD_2022.tse)
#' 
getFeatureStats <- function(tse, assay.type = NULL){
  
  if(is.null(assay.type)){
    assay.type <- assayNames(tse)[1]
    message(sprintf("assay.type defaulting to: %s", assay.type))
  }
  
  assay.obj <- assay(tse, assay.type)
  
  tse_stats.df <- as.data.frame(t(apply(assay.obj, 1, function(x)
    return(c(
      N = length(x), 
      N_zero = sum(x == 0), 
      N_not_zero = sum(x != 0), 
      Mean = mean(x), 
      Median = median(x), 
      Min = min(x), 
      Max = max(x)
    )))))
  
  return(tse_stats.df)
}