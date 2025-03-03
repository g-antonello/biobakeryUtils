#' wrangle a metaphlan profiles table
#'
#' @param mpa a Metaphlan matrix as a data frame
#' @param taxonomic_lvl \code{character} indicating the taxonomic level wanted. Default is `Species`
#'
#' @return a \code{list containing both `profiles` and `taxonomies`}
#' @importFrom purrr transpose reduce
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#' 
#' data("metaphlanTestData")
#' wrangled_metaphlans <- sapply(c("Kingdom", "Class","Species"),
#' function(taxLvl)
#' wrangle_metaphlan(metaphlanTestData, taxonomic_lvl = taxLvl), 
#' USE.NAMES = TRUE, simplify = FALSE)

wrangle_metaphlan <- function(mpa, taxonomic_lvl = "Species"){
  
  taxonomies <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
  
  # for each row, find to which depth of taxonomy it arrives (as integers)
  # the +1 at the end ensures that kingdom is level 1 and not 0
  tax_lengths <- lengths(regmatches(mpa$clade_name, gregexpr("\\|", mpa$clade_name))) + 1
  tax_lvl_int <- which(tolower(taxonomic_lvl) == tolower(taxonomies))

  selection_vec <- tax_lengths == tax_lvl_int
  
  # initialize taxa that have been added
  taxa_added <- list()
  
  for(lvl in tax_lvl_int:2){
    # separate the taxa that have been chosen...
    mpa_lvlChosen <- mpa$clade_name[selection_vec]
    # from those that haven't
    mpa_leftovers <- mpa$clade_name[!selection_vec]
    
    # find out the taxonomic level at which the leftovers stop, this is analogous 
    # to what I did at the beginning of the function, but on the subset
    tax_lengths <- lengths(regmatches(mpa_leftovers, gregexpr("\\|", mpa_leftovers))) + 1
    
    # this is an important step: check in the 'leftovers' at the taxonomic level 
    # above the one of interest (eg. if species is the one of interest, it checks Genus).
    # then it does the same at the higher taxonomic level again
    taxa_check <- sapply(strsplit(mpa_leftovers[tax_lengths == (lvl - 1)], "\\|"),  "[", lvl - 1)
    # now check if among these taxa at one level higher, any was left behind in
    # the selection done at lower taxonomic level (i.e. 'which are the genera that
    # were not resolved at species level?')
    taxa_to_add <- sapply(taxa_check, function(tx) grep(tx, mpa_lvlChosen, value = TRUE)) |> 
      sapply(length)
    
    if(any(taxa_to_add == 0)){
      taxa_to_add <- taxa_to_add[taxa_to_add == 0]
      
      # update the selection vector with the newly added taxa
      taxa_to_add <- mpa_leftovers[grepl(names(taxa_to_add), mpa_leftovers)]
      
      selection_vec <- selection_vec | (mpa$clade_name == taxa_to_add)  
      # save iterations in a list that can be used later
      taxa_added[[taxonomies[lvl - 1]]] <- taxa_to_add
    }
    
  }
  
  
  mpa_refined <- as.data.frame(mpa[selection_vec,])
  
  # perform these actions only if there are actually taxa to add
  if(length(taxa_added) > 0){
    unknown_taxa_names_completed <-  sapply(
      Reduce(c, taxa_added),
      complete_unknown_taxonomy,
      tax_lvl_int = tax_lvl_int,
      simplify = F,
      USE.NAMES = T
    ) |>
      unlist()
    
    #' replace original names with those that were gradually added going up the
    #' taxonomic tree
    
    mpa_refined$clade_name <- sapply(mpa_refined$clade_name, function(x) 
      ifelse(
        x %in% names(unknown_taxa_names_completed),
        unknown_taxa_names_completed[x],
        x
      )
    )
  }
  
  taxonomy.df <- separate(
    select(mpa_refined["clade_name"], clade_name),
    col = clade_name,
    sep = "\\|",
    into = taxonomies[1:tax_lvl_int],
    fill = "right"
    )
  rownames(taxonomy.df) <- taxonomy.df[[tax_lvl_int]]
  
  # finish formatting mpa refined
  rownames(mpa_refined) <- sapply(strsplit(mpa_refined$clade_name, "\\|"), "[", tax_lvl_int)
  mpa_refined$clade_name <- NULL
  
  
  return(
    list(
      profiles = mpa_refined,
      taxonomies = taxonomy.df
    )
  )
}
