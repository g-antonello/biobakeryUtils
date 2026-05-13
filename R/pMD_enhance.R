###############################################################################
# some functions to format parkinsonsMetagenomicData data that users get from 
# loadParquetData or returnSamples
###############################################################################

#' Enhance parkinsonsMetagenomicData output
#' 
#' Formats and adds pieces to select outputs of parkinsonsMetagenomicData
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of
#'  \code{loadMetagenomcData} or \code{returnSamples}
#' @param data_type \code{character}. `data_type` parameter used in 
#'  \code{returnSamples}
#' 
#' @returns \code{TreeSummarizedExperiment}
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay colData colData<-
#' @export 
#'
#' @examples
#' 
#' suppressMessages(library(TreeSummarizedExperiment))
#' data_types <- c("relative_abundance",
#' "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' # MetaPhlAn
#' data(Bedarf_pMD_raw_MetaPhlAn.tse)
#' 
#' tse_enhanced <- pMD_enhance(Bedarf_pMD_raw_MetaPhlAn.tse, 
#'   data_type = data_types[1])
#'   
#' # HUMAnN pathways
#' data(Bedarf_pMD_raw_HUMAnN_pwy.tse)
#' 
#' tse_enhanced <- pMD_enhance(Bedarf_pMD_raw_HUMAnN_pwy.tse, 
#'   data_type = data_types[2])
#'   
#' # HUMAnN Gene Families - example takes too long, see example in 
#' # `pMD_enhance_HUMAnN_genefam`
#' # data(Bedarf_pMD_raw_HUMAnN_GeneFam.tse)
#' 
#' # tse_enhanced <- pMD_enhance(Bedarf_pMD_raw_HUMAnN_GeneFam.tse, 
#' #   data_type = data_types[3])
pMD_enhance <- function(input.tse, data_type){
  if(!any(data_type %in% c("relative_abundance", "pathabundance_unstratified", "genefamilies_unstratified"))){
    warning(sprintf("%s is not yet supported, returning unmodified input", data_type))
    return(input.tse)
  } else{
    
    # call specific function
    if(data_type == "relative_abundance"){
      return(pMD_enhance_MetaPhlAn(input.tse))
    }
    
    if(data_type == "pathabundance_unstratified"){
      return(pMD_enhance_HUMAnN_pwy(input.tse))
    }
    
    if(data_type == "genefamilies_unstratified"){
      return(pMD_enhance_HUMAnN_genefam(input.tse))
    }
  }
}

###############################################################################

#' Enhance MetaPhlAn data
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of 
#' \code{loadMetagenomcData} or \code{returnSamples}. works exclusively for the
#' \code{data_type = "relative_abundance"} case
#' 
#' @returns \code{TreeSummarizedExperiment}
#' 
#' @importFrom SummarizedExperiment rowData colData assay assays assayNames
#' @importFrom SummarizedExperiment rowData<- colData<- assay<- assays<- 
#' @importFrom SummarizedExperiment assayNames<- 
#' @importFrom S4Vectors metadata metadata<-
#' 
#' @export
#'
#' @examples
#' suppressMessages(library(TreeSummarizedExperiment))
#' data_types <- c("relative_abundance",
#' "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' data(Bedarf_pMD_raw_MetaPhlAn.tse)
#' 
#' tse_enhanced <- pMD_enhance_MetaPhlAn(Bedarf_pMD_raw_MetaPhlAn.tse)
#' 
#' tse_enhanced

pMD_enhance_MetaPhlAn <- function(input.tse){
  
  # fix potential NAs that arise when parkinsonsMetagenomicData retrieves
  # data and merges assays together
  assay(input.tse)[is.na(assay(input.tse))] <- 0
  
  
  # Fix cases of columns that contain tab character, replace '\t' with '  '
  # as far as I know, only metaphlan headers have this issue
  colData(input.tse) <- DataFrame(apply(colData(input.tse), 2, function(x) gsub("\t", "  ", x)))
  
  # filter input.tse to have only t__SGB level data, all columns of the Assay should sum to 100 at the 5th precision digit
  # first add the exception of UNCLASSIFIED, which should be added as it is to the lowest levels
  input.tse <- input.tse[!is.na(rowData(input.tse)$clade_name_terminal) | rowData(input.tse)$clade_name_kingdom == "UNCLASSIFIED",]
  
  # then keep only the first useful columns and rename them, put the rest in metadata
  new_rowData <- as.data.frame(rowData(input.tse))
  
  new_rowData_wanted_part <- new_rowData[,which(colnames(new_rowData) == "clade_name_kingdom"):which(colnames(new_rowData) == "clade_name_terminal")]
  colnames(new_rowData_wanted_part) <- all_taxonomy_levels  
  
  # Fill UNCLASSIFIED row all with UNCLASSIFIED
  new_rowData_wanted_part["UNCLASSIFIED",] <- paste(c("k", "p", "c", "o", "f", "g", "s", "t"), "UNCLASSIFIED", sep = "__")
  
  # add new rowData to input.tse and make shorter rownames
  rowData(input.tse) <- DataFrame(new_rowData_wanted_part)
  rownames(input.tse) <- rowData(input.tse)[["SGB"]]
  
  # create the extra_rowData data.frame to go into the metadata
  extra_rowData <- new_rowData[, which(grepl("NCBI", colnames(new_rowData)))]
  rownames(extra_rowData) <- rownames(input.tse)
  metadata(input.tse)$rowData_extraColumns <- extra_rowData
  
  # do some cleaning of useful metaphlan run data
  input.tse@colData$reads_processed <- as.integer(gsub("#([0-9]+).*", "\\1", input.tse@colData$reads_processed))
  input.tse@colData$db_version <- gsub("#", "", input.tse@colData$db_version, fixed = TRUE)
  
  # clean metaphlan command function
  colnames(input.tse@colData)[which(grepl("command", colnames(input.tse@colData)))] <- "metaphlan_command"
  
  # calculate and organize alternative useful assays
  assayNames(input.tse) <- "percent"
  assay(input.tse, "relative_abundance") <- assay(input.tse)/100
  assay(input.tse, "counts") <- sweep(assay(input.tse, "relative_abundance"),2, input.tse@colData$reads_processed, FUN =  "*")
  
  # reorder the assays so that relative abundance [0.1]
  # comes before the percent (default MetaPhlAn output)
  assays(input.tse) <- assays(input.tse)[c(2,1,3)]
  
  return(input.tse)
}

###############################################################################

#' Enhance HUMAnN unstratified pathways data
#'
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of 
#' \code{loadMetagenomcData} or \code{returnSamples}. works exclusively for the
#' \code{data_type = "pathabundance_unstratified"} case
#' 
#' @returns \code{TreeSummarizedExperiment}
#' 
#' @importFrom SummarizedExperiment rowData colData assay assays assayNames
#' @importFrom SummarizedExperiment rowData<- colData<- assay<- assays<- 
#' @importFrom SummarizedExperiment assayNames<- 
#' @importFrom S4Vectors metadata metadata<-
#' 
#' @export
#'
#' @examples
#' suppressMessages(library(TreeSummarizedExperiment))
#' data_types <- c("relative_abundance",
#' "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' data(Bedarf_pMD_raw_HUMAnN_pwy.tse)
#' 
#' tse_enhanced <- pMD_enhance_HUMAnN_pwy(Bedarf_pMD_raw_HUMAnN_pwy.tse)
#' 
#' tse_enhanced

pMD_enhance_HUMAnN_pwy <-function(input.tse){
  
  # fix potential NAs that arise when parkinsonsMetagenomicData retrieves
  # data and merges assays together
  assay(input.tse)[is.na(assay(input.tse))] <- 0
  
  # Fix cases of columns that contain tab character, replace '\t' with '  '
  # as far as I know, only metaphlan headers have this issue
  colData(input.tse) <- DataFrame(apply(colData(input.tse), 2, function(x) gsub("\t", "  ", x)))
  
  # fix rowData and rownames a little
  new_rowData <- cbind.data.frame(rowData(input.tse)$pathway, do.call(rbind, strsplit(x = rowData(input.tse)$pathway, ": ", fixed = TRUE)))
  colnames(new_rowData) <- c("Pathway", "MetaCyc_code", "Description")
  new_rowData$MetaCyc_code_safeName <- make.names(new_rowData[,"MetaCyc_code"])
  new_rowData <- new_rowData[,c(1,2,4,3)]
  rownames(new_rowData) <- new_rowData[["Pathway"]]
  
  # add new rowData
  rowData(input.tse) <- DataFrame(new_rowData)
  # rename rownames for better compatibility
  rownames(input.tse) <- rowData(input.tse)$MetaCyc_code_safeName
  
  # re-derive transformations
  assay(input.tse, "relative_abundance") <- apply(assay(input.tse), 2, function(x) x/sum(x))
  assay(input.tse, "cpm") <- apply(assay(input.tse, "relative_abundance"), 2, function(x) x*10^6)
  
  return(input.tse)
}

################################################################################

#' Enhance HUMAnN unstratified gene families data
#' 
#' @param input.tse \code{TreeSummarizedExperiment}. as it comes out of 
#' \code{loadMetagenomcData} or \code{returnSamples}. works exclusively for the
#' \code{data_type = "genefamilies_unstratified"} case
#' 
#' @returns \code{TreeSummarizedExperiment}
#' 
#' @importFrom SummarizedExperiment rowData colData assay assays assayNames
#' @importFrom SummarizedExperiment rowData<- colData<- assay<- assays<- 
#' @importFrom SummarizedExperiment assayNames<- 
#' 
#' @export
#'
#' @examples
#' # this example takes a long time, but it's the problem of Gene Families:
#' # They are too many.
#' 
#' # library(parkinsonsMetagenomicData)
#' # data("sampleMetadata", package = "parkinsonsMetagenomicData")
#' 
#' # data_types <- c("relative_abundance",
#' #                "pathabundance_unstratified", "genefamilies_unstratified")
#' 
#' # example_metadata.df <- sampleMetadata[grepl("Bedarf",
#' #                                             sampleMetadata$study_name),]
#' 
#' # tmp.tse <- returnSamples(example_metadata.df, data_type = data_types[3])
#' 
#' # Bedarf_pMD_raw_HUMAnN_GeneFam.tse <- pMD_enhance_HUMAnN_genefam(tmp.tse)
#' 
#' # Bedarf_pMD_raw_HUMAnN_GeneFam.tse

pMD_enhance_HUMAnN_genefam <-function(input.tse){
  
  # fix potential NAs that arise when parkinsonsMetagenomicData retrieves
  # data and merges assays together
  assay(input.tse)[is.na(assay(input.tse))] <- 0
  
  # Fix cases of columns that contain tab character, replace '\t' with '  '
  # as far as I know, only metaphlan headers have this issue
  colData(input.tse) <- DataFrame(apply(colData(input.tse), 2, function(x) gsub("\t", "  ", x)))
  
  # re-derive transformations
  assay(input.tse, "relative_abundance") <- apply(assay(input.tse), 2, function(x) x/sum(x))
  assay(input.tse, "cpm") <- apply(assay(input.tse, "relative_abundance"), 2, function(x) x*10^6)
  
  return(input.tse)
}
