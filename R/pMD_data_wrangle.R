###############################################################################
# some functions to format parkinsonsMetagenomicData data that users get from 
# loadParquetData or returnSamples
###############################################################################

#' Enhance parkinsonsMetagenomicData output
#' 
#' Formats and adds pieces to select outputs of parkinsonsMetagenomicData
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of \code{loadMetagenomcData}
#' or \code{returnSamples}
#' @param sampleMetadata.df \code{data.frame}. A data.frame to swap the default 
#' with. NB: it must have a `uuid` column with valid uuids
#' \code{loadMetagenomcData} or \code{returnSamples}
#' @param data_type \code{character}. that got retrieved using 
#' \code{loadMetagenomcData} or \code{returnSamples}
#' @param addPhyloTree \code{logical}. Should a phylogenetic `rowTree` be added 
#' to the \code{TreeSummarizedExperiment}? Works only with MetaPhlAn input (`data_type = "relative_abundance"`)
#' 
#' @returns \code{TreeSummarizedExperiment}
#' @export 
#'
#' @examples
#' \dontrun{
#' library(parkinsonsMetagenomicData)
#' 
#' data_types <- c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' hf_con <- accessParquetData(data_types = data_types)
#' 
#' example_metadata.df <- dplyr::filter(sampleMetadata, grepl("Bedarf", study_name))
#' 
#' 
#' # MetaPhlAn
#' 
#' tse_basic <- loadParquetData(con = hf_con, data_type = "relative_abundance", filter_values = list(uuid = example_metadata.df$uuid))
#' 
#' tse_enhanced <- pMD_enhance(tse_basic, sampleMetadata.df = example_metadata.df, data_type = "relative_abundance", addPhyloTree = TRUE)
#' 
#' # HUMAnN pathways
#' tse_basic <- loadParquetData(con = hf_con, data_type = "pathabundance_unstratified", filter_values = list(uuid = example_metadata.df$uuid))
#' 
#' tse_enhanced <- pMD_enhance(tse_basic, sampleMetadata.df = example_metadata.df, data_type = "pathabundance_unstratified", addPhyloTree = TRUE)
#' 
#' # HUMAnN gene families
#' tse_basic <- loadParquetData(con = hf_con, data_type = "genefamilies_unstratified", filter_values = list(uuid = example_metadata.df$uuid))
#' 
#' tse_enhanced <- pMD_enhance(tse_basic, sampleMetadata.df = example_metadata.df, data_type = "genefamilies_unstratified", addPhyloTree = TRUE)
#' 
#' }

pMD_enhance <- function(input.tse, sampleMetadata.df, data_type, addPhyloTree = FALSE){
  if(!any(data_type %in% c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified"))){
    warning(sprintf("%s is not yet supported, returning unmodified input", data_type))
    return(input.tse)
  } else{
    
    # fix potential NAs that arise when parkinsonsMetagenomicData retrieves
    # data and merges assays together
    assay(input.tse)[is.na(assay(input.tse))] <- 0
    
    # call specific function
    if(data_type == "relative_abundance"){
      return(pMD_enhance_MetaPhlAn(input.tse, sampleMetadata.df, addPhyloTree = addPhyloTree))
    }
    
    if(data_type == "pathabundance_unstratified"){
      return(pMD_enhance_HUMAnN_pwy(input.tse, sampleMetadata.df))
    }
    
    if(data_type == "genefamilies_unstratified"){
      return(pMD_enhance_HUMAnN_genefam(input.tse, sampleMetadata.df))
    }
  }
}

###############################################################################

#' Enhance MetaPhlAn data
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of \code{loadMetagenomcData}
#' or \code{returnSamples} and `data_type = 'relabundance'`
#' @param sampleMetadata.df \code{data.frame}. A data.frame to swap the default 
#' with. NB: it must have a `uuid` column with valid uuids
#' @param addPhyloTree \code{logical}. Should a phylogenetic `rowTree` be added 
#' to the \code{TreeSummarizedExperiment}? Works only with MetaPhlAn input (`data_type = "relative_abundance"`)
#' 
#' @returns \code{TreeSummarizedExperiment}
#' @export
#'
#' @examples
#' \dontrun{
#' library(parkinsonsMetagenomicData)
#' 
#' data_types <- c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' hf_con <- accessParquetData(data_types = data_types)
#' 
#' example_metadata.df <- dplyr::filter(sampleMetadata, grepl("Bedarf", study_name))
#' 
#' 
#' # MetaPhlAn
#' 
#' tse_basic <- loadParquetData(con = hf_con, data_type = "relative_abundance", filter_values = list(uuid = example_metadata.df$uuid))
#' 
#' tse_enhanced <- pMD_enhance(tse_basic, sampleMetadata.df = example_metadata.df, data_type = "relative_abundance", addPhyloTree = TRUE)
#' }

pMD_enhance_MetaPhlAn <- function(input.tse, sampleMetadata.df, addPhyloTree = TRUE){
  
  # filter input.tse to have only t__SGB level data, all columns of the Assay should sum to 100 at the 5th precision digit
  # first add the exception of UNCLASSIFIED, which should be added as it is to the lowest levels
  input.tse <- input.tse[!is.na(rowData(input.tse)$clade_name_terminal) | rowData(input.tse)$clade_name_kingdom == "UNCLASSIFIED",]
  
  # then keep only the first useful columns and rename them, put the rest in metadata
  new_rowData <- rowData(input.tse) %>% 
    as.data.frame()
  
  new_rowData_wanted_part <- select(new_rowData, clade_name_kingdom:clade_name_terminal) %>% 
    dplyr::rename_all(.funs = function(x) {
      gsub("clade_name_", "", x) %>% 
        tools::toTitleCase()
    }) %>% 
    dplyr::rename(SGB = Terminal)
  
  # Fill UNCLASSIFIED row all with UNCLASSIFIED
  new_rowData_wanted_part["UNCLASSIFIED",] <- paste(c("k", "p", "c", "o", "f", "g", "s", "t"), "UNCLASSIFIED", sep = "__")
  
  # add new rowData to input.tse and make shorter rownames
  rowData(input.tse) <- DataFrame(new_rowData_wanted_part)
  rownames(input.tse) <- rowData(input.tse)[["SGB"]]
  
  # create the extra_rowData data.frame to go into the metadata
  extra_rowData <- select(new_rowData, !(clade_name_kingdom:clade_name_terminal))
  rownames(extra_rowData) <- rownames(input.tse)
  metadata(input.tse)$rowData_extraColumns <- extra_rowData
  
  # do some cleaning of useful metaphlan run data
  input.tse@colData$number_reads <- parse_number(input.tse@colData$reads_processed)
  input.tse@colData$reads_processed <- NULL
  input.tse@colData$db_version <- gsub("#", "", input.tse@colData$db_version, fixed = TRUE)
  
  # clean metaphlan command function
  colnames(input.tse@colData)[which(grepl("command", colnames(input.tse@colData)))] <- "metaphlan_command"
  
  # calculate and organize alternative useful assays
  assayNames(input.tse) <- "percent"
  assay(input.tse, "relabundance") <- assay(input.tse)/100
  assay(input.tse, "counts") <- t(apply(assay(input.tse, "relabundance"), 1, function(x) x * input.tse@colData$number_reads))
  
  # reorder the assays so that relative abundance [0.1]
  # comes before the percent (default MetaPhlAn output)
  assays(input.tse) <- assays(input.tse)[c(2,1,3)]
  if(addPhyloTree){
    cphl_version <- unique(input.tse@colData$db_version)
    cphl_version <- as.character(parse_number(gsub(".*PhlAnSGB_", "", cphl_version)))
    
    if(length(cphl_version) == 1){
      # add phylogenetic tree
      input.tse <-  AddPhyloTree_to_mpa_tse(data.tse = input.tse, CHOCOPhlAn_version = cphl_version)  
    } else{
      message("cannot assign a tree to features that are not called under the same CHOCOPhlAn version")
    }
  }
  
  return(input.tse)
}

###############################################################################

#' Enhance HUMAnN unstratified pathways data
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of 
#' \code{loadMetagenomcData} or \code{returnSamples}. using 
#' `data_type = "pathabundance_unstratified"`
#' @param sampleMetadata.df \code{data.frame}. A data.frame to swap the default 
#' with. NB: it must have a `uuid` column with valid uuids
#' 
#' @returns \code{TreeSummarizedExperiment}
#' @export
#'
#' @examples
#' \dontrun{
#' library(parkinsonsMetagenomicData)
#' 
#' data_types <- c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' hf_con <- accessParquetData(data_types = data_types)
#' 
#' example_metadata.df <- dplyr::filter(sampleMetadata, grepl("Bedarf", study_name))
#' 
#' # HUMAnN pathways
#' tse_basic <- loadParquetData(con = hf_con, data_type = "pathabundance_unstratified", filter_values = list(uuid = example_metadata.df$uuid))
#' 
#' tse_enhanced <- pMD_enhance(tse_basic, sampleMetadata.df = example_metadata.df, data_type = "pathabundance_unstratified", addPhyloTree = TRUE)
#' 
#' }

pMD_enhance_HUMAnN_pwy <-function(input.tse, sampleMetadata.df){
  
  # fix rowData and rownames a little
  new_rowData <- rowData(input.tse) %>% 
    as.data.frame() %>%
    separate_wider_delim(cols = "pathway", delim = ": ", names = c("MetaCyc_code", "Description"), cols_remove = FALSE, too_few = "align_start") %>% 
    mutate(MetaCyc_code_safeName = make.names(MetaCyc_code)) %>% 
    relocate(pathway)
  rowData(input.tse) <- DataFrame(new_rowData)
  # rename rownames for better compatibility
  rownames(input.tse) <- rowData(input.tse)$MetaCyc_code_safeName
  
  # re-derive transformations
  assay(input.tse, "relabundance") <- apply(assay(input.tse), 2, function(x) x/sum(x))
  assay(input.tse, "cpm") <- apply(assay(input.tse, "relabundance"), 2, function(x) x*10^6)
  
  return(input.tse)
}

################################################################################

#' Enhance HUMAnN unstratified gene families data
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of 
#' \code{loadMetagenomcData} or \code{returnSamples}. using 
#' `data_type = "genefamilies_unstratified"`
#' @param sampleMetadata.df \code{data.frame}. A data.frame to swap the default 
#' with. NB: it must have a `uuid` column with valid uuids
#' 
#' @returns \code{TreeSummarizedExperiment}
#' @export
#'
#' @examples
#' \dontrun{
#' library(parkinsonsMetagenomicData)
#' 
#' data_types <- c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' hf_con <- accessParquetData(data_types = data_types)
#' 
#' example_metadata.df <- dplyr::filter(sampleMetadata, grepl("Bedarf", study_name))
#' 
#' # HUMAnN gene families
#' tse_basic <- loadParquetData(con = hf_con, data_type = "genefamilies_unstratified", filter_values = list(uuid = example_metadata.df$uuid))
#' 
#' tse_enhanced <- pMD_enhance(tse_basic, sampleMetadata.df = example_metadata.df, data_type = "genefamilies_unstratified", addPhyloTree = TRUE)
#' 
#' }

pMD_enhance_HUMAnN_genefam <-function(input.tse, sampleMetadata.df){
  
  # re-derive transformations
  assay(input.tse, "relabundance") <- apply(assay(input.tse), 2, function(x) x/sum(x))
  assay(input.tse, "cpm") <- apply(assay(input.tse, "relabundance"), 2, function(x) x*10^6)
  
  return(input.tse)
}

#' Replace colData of a (Tree)SummarizedExperiment
#' 
#' A pipe-friendly wrapper to replace colData from a (Tree)SummarizedExperiment.
#' By design, if a new colData with more samples is provided, only those that map
#' to input.se are inserted in the colData slot.
#'
#' @param input.se \code{(Tree)SummarizedExperiment} 
#' @param new.df \code{data.frame or DataFrame} with rownames corresponding to 
#' rownames in `rownames(colData(input.se))`. The correct order is coerced.
#'
#' @returns \code{(Tree)SummarizedExperiment} with replaced colData
#' @export
#'
#' @examples
#' \dontrun{
#' library(parkinsonsMetagenomicData)
#' 
#' data_types <- c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' hf_con <- accessParquetData(data_types = data_types)
#' 
#' example_metadata.df <- dplyr::filter(sampleMetadata, grepl("Bedarf", study_name))
#' 
#' example_metadata.df$ABC_Var <- factor(sample(1:0, nrow(example_metadata.df), replace = TRUE))
#' rownames(example_metadata.df) <- example_metadata.df$uuid
#' 
#' # HUMAnN gene families
#' tse_basic <- loadParquetData(
#'   con = hf_con, 
#'   data_type = "genefamilies_unstratified", 
#'   filter_values = list(uuid = example_metadata.df$uuid)
#' )
#' 
#' "ABC_Var" %in% colnames(colData(tse_basic))
#' 
#' tse_newColData <- tse_basic %>%
#'   replace_colData(  new.df = example_metadata.df)
#'
#' "ABC_Var" %in% colnames(colData(tse_newColData))
#'   
#' }

replace_colData <- function(input.se, new.df){
  
  if(!all(rownames(colData(input.se)) %in% rownames(new.df))){
    stop("some rownames of old colData are missing from the new colData")
  }
  
  colData(input.se) <- DataFrame(new.df[rownames(colData(input.se)),])
  return(input.se)
}
