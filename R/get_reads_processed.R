

#' Get reads processed by MetaPhlAn 3 and 4
#' 
#' A simple function to read in the .gz files cached with 
#' `parkinsonsMetagenomicData::cacheMetagenomicData()`
#' 
#' @param metaphlan_profile_paths \code{character}, The path to the metaphlan 
#' profile.
#' 
#' @returns A \code{integer} vector containing the number of reads processed
#' 
#' @importFrom  readr parse_number
#' 
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' # this doesn't work yet, at some point I will put an example file into 
#' # extdata
#'   path_to_mpa <- system.file("mpa_single_file.tsv.gz", "extdata", 
#'   package = "biobakeryUtils")
#'   
#'   get_reads_processed(path_to_mpa)
#' }

get_reads_processed <- function(metaphlan_profile_paths){
  # make a length-1 vector function
  
  get_reads_processed.single <- function(metaphlan_profile_path){
    return(parse_number(readLines(metaphlan_profile_path, n = 3)[[3]]))
  }
  
  return(as.integer(sapply(
    metaphlan_profile_paths, get_reads_processed.single
  )))
}
