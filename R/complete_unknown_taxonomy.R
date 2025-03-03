#' Complete Unknown taxonomy
#'  
#' This is a utility function for the `metaphlan_wrangle function`, which needs
#' it in case it finds taxa that are unique at higher taxonomies levels than the 
#' one requested (e.g. a Genus without any species mapped to it. maybe it is still 
#' a valuable genus, so it has to be repeated with the same name at the taxonomic level below).
#' This is also valuable because another approach would discard it, and the total
#' relative abundances estimated would not add up.
#' 
#' @param x \code{character},  the taxonomy as it appears in the first column of a metaphlan table
#' @param tax_lvl_int \code{character} the level at which the user wants the taxonomies aggregated. The default goes all the way to the SGB level (`8`). Species level is `7`, Genus is `6` and so on.
#'
#' @return \code{character}, the same as the input vector, with the lowest taxonomy repeated down to the `tax_lvl_int` chosen
#' @importFrom dplyr separate_wider_delim
#' @export
#'
#' @examples
#' 
#' example_taxonomies <- c("k__Bacteria|p__Firmicutes|c__CFGB1347", "k__Bacteria|p__Actinobacteria", "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Eubacteriales_Family_XIII_Incertae_Sedis|g__Lentihominibacter|s__Lentihominibacter_faecis")
#' 
#' # complete taxonomies to species level (level of taxonomic depth)
#' complete_unknown_taxonomy(example_taxonomies, tax_lvl_int = 7)
#' 
#' 

complete_unknown_taxonomy <- function(x, tax_lvl_int = 8){
  
  #input checks 
  if (!is.character(x)) {
    stop("Input 'x' must be a character vector.")
  }
  
  if (!is.numeric(tax_lvl_int) || tax_lvl_int < 1 || tax_lvl_int > 8 || tax_lvl_int != round(tax_lvl_int)) {
    stop("'tax_lvl_int' must be an integer between 1 and 8.")
  }
  
  tmp <- sapply(strsplit(x, "\\|"), function(l) l[length(l)])
  possible_levels <- c("k__","p__","c__","o__", "f__", "g__", "s__", "SGB_")[1:tax_lvl_int]
  latest_taxonomy <- which(grepl(substr(tmp, 1, 3), possible_levels))
  if(length(latest_taxonomy) == 0){
    latest_taxonomy <- 0
    return(paste(rep(tmp, tax_lvl_int), collapse = "|"))
  } else{
    return(paste(c(x, rep(tmp, tax_lvl_int - latest_taxonomy)), collapse = "|")) 
  }
}