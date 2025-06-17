#' A vector of taxonomies and their typical prefixes in metaphlan
#'
#' This vector contains taxonomies from Kingdom to SGB level. the vector is
#' generated as : all_taxonomy_levels <- c( k__ = "Kingdom", p__ = "Phylum",
#'   c__ = "Class", o__ = "Order", f__ = "Family", g__ = "Genus", 
#'  s__ = "Species", t__ = "SGB")
# 
# save(all_taxonomy_levels, file = "data/all_taxonomy_levels.rda")

#' 
#' @format A character vector with 8 elements:
#' \describe{
#'   \item{all_taxonomy_levels}{The taxonomy level with capitalized 1st letter}
#' }
"all_taxonomy_levels"

