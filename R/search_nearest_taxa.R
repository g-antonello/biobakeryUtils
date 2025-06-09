
#' Find phylogenetically closest taxa to a target
#'
#' @param phyloTree \code{phylo} tree that can be handled by `ape` and 
#' `TreeTools` CRAN packages. IMPORTANT: So far, this function was thought to
#'  work on full taxonomy as tip labels. So before extracting the phyloTree,
#'  please run `rename_rownames.tse(your_data.tse, collapse_taxonomy(as.data.frame(rowData(your_data.tse))))`
#' @param target_SGB \code{character} with SGB name either as "t__SGB12345" or 
#' as full taxonomy.
#' @param taxonomy_sep \code{character} specifying the separation character used
#' in the tree's tip labels. Default is "|"
#'
#' @importFrom ape cophenetic.phylo
#' @returns
#' @export
#'
#' @examples

search_nearest_taxa <- function(phyloTree, target_SGB, taxonomy_sep = "|"){

  # get pariwise distances
  ape_coph_distances <- cophenetic.phylo(phyloTree)
  
  # remove the self-comparison
  diag(ape_coph_distances) <- NA
  
  # get only the comparisons against the target taxon
  all_dists_vs_target <- sort(ape_coph_distances[grepl(target_SGB, rownames(ape_coph_distances)),])
  
  intermediate.df <- expand_taxonomy(names(all_dists_vs_target), sep = taxonomy_sep) %>% 
  mutate(
    phylo_dist = all_dists_vs_target,
    known_species = !grepl("GGB|SGB|_bacterium|UNCLASSIFIED", Species) | grepl("Isolate|isolate", Species)
  ) %>% 
    arrange(phylo_dist)
  
  return(intermediate.df)
}
