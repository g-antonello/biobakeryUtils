#' Complete Unknown taxonomy
#'  
#' This is a utility function for the `metaphlan_wrangle function`, which needs
#' it in case it finds taxa that are unique at higher taxonomies levels than the 
#' one requested (e.g. a Genus without any species mapped to it. maybe it is 
#' still a valuable genus, so it has to be repeated with the same name at the 
#' taxonomic level below). This is also valuable because another approach would
#' discard it, and the total relative abundances estimated would not add up.
#' 
#' @param x \code{character},  the taxonomy as it appears in the first column of
#' a metaphlan table
#' @param tax_lvl_int \code{character} the level at which the user wants the 
#' taxonomies aggregated. The default goes all the way to the SGB level (`8`).
#' Species level is `7`, Genus is `6` and so on.
#' @return \code{character}, the same as the input vector, with the lowest taxonomy repeated down to the `tax_lvl_int` chosen
#' @importFrom tidyr separate_wider_delim
#' @importFrom dplyr select
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
  possible_levels <- c("k__","p__","c__","o__", "f__", "g__", "s__", "t__")[1:tax_lvl_int]
  latest_taxonomy <- which(grepl(substr(tmp, 1, 3), possible_levels))
  if(length(latest_taxonomy) == 0){
    latest_taxonomy <- 0
    return(paste(rep(tmp, tax_lvl_int), collapse = "|"))
  } else{
    return(paste(c(x, rep(tmp, tax_lvl_int - latest_taxonomy)), collapse = "|")) 
  }
}


#' wrangle a metaphlan profiles table
#' 
#' Works really similarly to mia::importMetaphlan, the difference being that 
#' this function starts from data already loaded in R. Also, it does not 
#' handel ncbi_tax_id so well
#' 
#' @param mpa a Metaphlan matrix as a data frame as comes out after 
#' `merge_metaphlan.py`
#' @param taxonomic_lvl \code{character} indicating the taxonomic level wanted. 
#' Default is `Species`
#' @param row.names \code{character} of length 1. Either `'short'`, `'long'`, or 
#' `'none'`. Default is `'short'`, which returns the highest taxonomic 
#' resolution the input without taxonomies above. This is passed to 
#' `expand_taxonomy()`
#' @param sep \code{character} Passed to `expand_taxonomy()` to specift the 
#' taxonomy separator in the input character vector. Default is "|"
#' @return a \code{list containing both `profiles` and `taxonomies`}
#' @importFrom purrr transpose reduce
#' @importFrom tidyr separate_wider_delim 
#' @importFrom dplyr select
#' @export
#'
#' @examples
#' 
#' wallen_profiles <- data.table::fread(system.file("extdata",
#'                                                  "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                                                  package = "biobakeryUtils"))
#' 
#' 
#' profiles <- sapply(c("Phylum", "Genus", "Species"), function(taxLvl) {
#'   wrangle_metaphlan(wallen_profiles, taxonomic_lvl = taxLvl)
#' },
#' USE.NAMES = TRUE, simplify = FALSE)
#' 
#' profiles$Species$profiles[1:5,1:5]
#' profiles$Species$taxonomies[1:5,]
#' 

wrangle_metaphlan <- function(mpa, taxonomic_lvl = "Species", row.names = "short", sep = "|"){
  
  # the package imports all_taxonomy_levels, a named character vector
  
  # for each row, find to which depth of taxonomy it arrives (as integers)
  # the +1 at the end ensures that kingdom is level 1 and not 0
  tax_lengths <- lengths(regmatches(mpa$clade_name, gregexpr("\\|", mpa$clade_name))) + 1
  tax_lvl_int <- which(tolower(taxonomic_lvl) == tolower(all_taxonomy_levels))

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
      taxa_added[[all_taxonomy_levels[lvl - 1]]] <- taxa_to_add
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
    
    # replace original names with those that were gradually added going up the
    # taxonomic tree
    
    mpa_refined$clade_name <- sapply(mpa_refined$clade_name, function(x) 
      ifelse(
        x %in% names(unknown_taxa_names_completed),
        unknown_taxa_names_completed[x],
        x
      )
    )
  }
  
  taxonomy.df <- expand_taxonomy(mpa_refined$clade_name, row.names = row.names, sep = sep)
  
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
