#' Specify classes of column data
#' Most importantly, the levels of a factor
#' 
#' @param df \code{data.frame} or \code{DataFrame} Object
#' 
#' @importFrom SummarizedExperiment colData
#' 
#' @returns A \code{data.frame} with per-column class and if a factor the levels 
#' in a third column
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#' 
#' # basic assay
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#' 
#' # sample data
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5), 
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#' 
#' # create rowData
#' rowData <- DataFrame(
#' Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#' Feature_gr2 = rep(c(2,3,4,6,1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#' 
#' # base tse
#' tse <- TreeSummarizedExperiment(
#'   assays = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#' 
#' # add altExp
#' alt_colData <- sampC
#' alt_colData$altTreat <- rep(0:1, 5)
#' 
#' alt_counts <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(assays = list(counts = alt_counts), rowData = alt_rowData, colData = alt_colData)
#' 
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#'  
#' DataFrameSpecs(colData(tse))
#' DataFrameSpecs(rowData(tse))

DataFrameSpecs <- function(df) {
  colSpecs.df <- data.frame(
    colName  = c("rownames", colnames(df)),
    colClass = c("character", vapply(df, function(x) class(x)[1], character(1))),
    fctLevels = c(NA_character_, vapply(df, function(x) {
      if (is.factor(x)) paste(levels(x), collapse = ";") else NA_character_
    }, character(1))),
    stringsAsFactors = FALSE
  )
  
  # readr specs
  colSpecs.df$readrClass <- substr(colSpecs.df$colClass, 1, 1)
  
  return(colSpecs.df)
}



#' Split a `TreeSummarizedExperiment` object into multiple human readable files
#'
#' The components are saved into a directory that gets created by this function.
#' AltExp slot is ignored in this function
#'
#' @param tse \code{TreeSummarizedExperiment} object with or without
#' rowTree and/or colTree
#' @param out.dir \code{character} specifying the directory in which a new directory
#' will be generated, with the name of the R object passed to the `tse`
#' parameter.
#' 
#' @import TreeSummarizedExperiment
#' @importFrom readr write_tsv
#' @importFrom jsonlite write_json
#' @importFrom ape write.tree
#' @importFrom S4Vectors metadata
#'
#' @export
#' 
#' @returns \code{NULL} - This function writes files on disk
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#' 
#' # basic assay
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#' 
#' # sample data
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5), 
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#' 
#' # create rowData
#' rowData <- DataFrame(
#' Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#' Feature_gr2 = rep(c(2,3,4,6,1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#' 
#' # base tse
#' tse <- TreeSummarizedExperiment(
#'   assays = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#' 
#' # add altExp
#' alt_colData <- sampC
#' alt_colData$altTreat <- rep(0:1, 5)
#' 
#' alt_counts <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(assays = list(counts = alt_counts), rowData = alt_rowData, colData = alt_colData)
#' 
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#' 
#' tmpdir <- file.path(tempdir(), "testTSE")
#' dir.create(tmpdir)
#' write_TSE_to_dir_noAltExp(tse, tmpdir)
#' list.files(tmpdir, recursive = TRUE)
#' unlink(tmpdir, recursive = TRUE)

write_TSE_to_dir_noAltExp <- function(tse, out.dir) {
  # write the colData
  write_tsv(
    colData(tse) |> as.data.frame() |> rownames_to_column("rownames"),
    file = file.path(out.dir, paste0("colData", ".tsv"))
  )
  # also write column specifications
  write_tsv(
    DataFrameSpecs(colData(tse)), file.path(out.dir, paste0("colData_colSpecs", ".tsv"))
  )
  
  # write the rowData
  write_tsv(
    rowData(tse) |> as.data.frame() |> rownames_to_column("rownames"),
    file = file.path(out.dir, paste0("rowData", ".tsv"))
  )
  # also write row specifications
  write_tsv(
    DataFrameSpecs(rowData(tse)), file.path(out.dir, paste0("rowData_colSpecs", ".tsv"))
  )
  
  # write metadata as json file
  jsonlite::write_json(metadata(tse), file.path(out.dir, "metadata.json"))
  
  # write the assays, including the sub-directory
  dir.create(
    file.path(out.dir, "assays"),
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  for (assayX in assayNames(tse)) {
    write_tsv(
      assay(tse, assayX) |>
        as.data.frame() |>
        rownames_to_column("rownames"),
      file = file.path(out.dir, "assays", paste0(assayX, ".tsv"))
    )
    
  }
  
  # write order of assays
  writeLines(assayNames(tse), con = file.path(out.dir, "assays", "assaysOrder.txt"), sep = "\n")
  
  # write rowTree file
  if (inherits(tse, "TreeSummarizedExperiment")) {
    if (!is.null(rowTree(tse))) {
      write.tree(rowTree(tse), file = file.path(out.dir, "rowTree.tre"))
    }
    # write colTree file
    if (!is.null(colTree(tse))) {
      write.tree(rowTree(tse), file = file.path(out.dir, "colTree.tre"))
    }
  }
}

#' Split a `TreeSummarizedExperiment` object into multiple human readable files
#'
#' The components are saved into a directory that gets created by this function.
#' AltExp slot is ignored in this function
#'
#' @param tse \code{TreeSummarizedExperiment} object with or without
#' rowTree and/or colTree
#' @param out.dir \code{character} specifying the directory in which a new directory
#' will be generated, with the name of the R object passed to the `tse`
#' parameter.
#'
#' @importFrom readr write_tsv
#' @importFrom jsonlite write_json
#' @importFrom tibble rownames_to_column
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom ape write.tree
#' @importFrom S4Vectors metadata
#'
#' @returns Nothing, the main output of this function is multiple .tsv files
#' @export
#' 
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#' 
#' # basic assay
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#' 
#' # sample data
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5), 
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#' 
#' # create rowData
#' rowData <- DataFrame(
#' Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#' Feature_gr2 = rep(c(2,3,4,6,1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#' 
#' # base tse
#' tse <- TreeSummarizedExperiment(
#'   assays = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#' 
#' # add altExp
#' alt_colData <- sampC
#' alt_colData$altTreat <- rep(0:1, 5)
#' 
#' alt_counts <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(assays = list(counts = alt_counts), rowData = alt_rowData, colData = alt_colData)
#' 
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#' 
#' tmpdir <- file.path(tempdir(), "testTSE_complete")
#' write_TSE_to_dir(tse, tmpdir)
#' list.files(tmpdir, recursive = TRUE)

write_TSE_to_dir <- function(tse, out.dir){
  if(dir.exists(out.dir)){
    stop("Output directory already exists")
  } 
  dir.create(out.dir, recursive = TRUE)
  
  write_TSE_to_dir_noAltExp(tse = tse, out.dir = out.dir)
  
  if(length(altExpNames(tse)) > 0){
  for(aexp in altExpNames(tse)){
    tse_aexp <- altExp(tse, aexp)
    
    altexp.outdir <- file.path(out.dir, "altExps", aexp)
    
    if(dir.exists(altexp.outdir)){
      stop("Output directory already exists")
    } 
    dir.create(altexp.outdir, recursive = TRUE)
    
    write_TSE_to_dir_noAltExp(tse_aexp, out.dir = altexp.outdir)
    }
    
  # also write loading order, possibly
  writeLines(altExpNames(tse), con = file.path(out.dir, "altExps", "AltExpOrder.txt"), sep = "\n")
  }
}

#' Read multiple tsv files into a TreeSummarizedExperiment
#'
#' This utility function is the counterpart of
#'
#' @param tse.dir \code{character} indicating the directory that was created
#' with `write_TSE_to_dir()`
#'
#' @returns A \code{TreeSummarizedExperiment}
#'
#' @importFrom data.table fread
#' @importFrom tibble column_to_rownames
#' @importFrom purrr reduce
#' @importFrom jsonlite read_json
#' 
#' @export
#' 
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#' 
#' # basic assay
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#' 
#' # sample data
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5), 
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#' 
#' # create rowData
#' rowData <- DataFrame(
#' Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#' Feature_gr2 = rep(c(2,3,4,6,1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#' 
#' # base tse
#' tse <- TreeSummarizedExperiment(
#'   assays = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#' 
#' # add altExp
#' alt_colData <- sampC
#' alt_colData$altTreat <- rep(0:1, 5)
#' 
#' alt_counts <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(assays = list(counts = alt_counts), rowData = alt_rowData, colData = alt_colData)
#' 
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#' 
#' tmpdir <- file.path(tempdir(), "testTSE3")
#' dir.create(tmpdir)
#' write_TSE_to_dir_noAltExp(tse, tmpdir)
#' list.files(tmpdir, recursive = TRUE)
#' 
#' # re-load tse
#' tse2 <- read_TSE_from_dir_noAltExp(tmpdir)
#' tse2
#' tse3 <- tse
#' altExp(tse3) <- NULL
#' all.equal(tse3, tse2)  # TRUE

read_TSE_from_dir_noAltExp <- function(tse.dir) {
  
  # colData
  ## get colClasses and specify them while reading colData
  colData_colSpecs <- fread(
    file.path(tse.dir, "colData_colSpecs.tsv"), sep = "\t", header = TRUE)
  
  colClasses <- colData_colSpecs$colClass
  names(colClasses) <- colData_colSpecs$colName
  
  ## read colData and specify colClasses
  colData <- fread(
    file = file.path(tse.dir, "colData.tsv"),
    sep = "\t",
    colClasses = colClasses
  ) |>
    column_to_rownames("rownames")
  
  ## recode colData factors with correct levels
  for(col in colData_colSpecs$colName[colData_colSpecs$colClass == "factor"]){
    colData[[col]] <- factor(colData[[col]], levels = strsplit(colData_colSpecs$fctLevels[colData_colSpecs$colName == col], ";|\\,\\ ")[[1]])
  }
  
  # rowData
  ## get colClasses and specify them while reading rowData
  rowData_colSpecs <- fread(
    file = file.path(tse.dir, "rowData_colSpecs.tsv"), sep = "\t", header = TRUE) 
  
  colClasses <- rowData_colSpecs$colClass
  names(colClasses) <- rowData_colSpecs$colName
  
  ## read rowData and specify colClasses
  rowData <- fread(
    file = file.path(tse.dir, "rowData.tsv"),
    colClasses = colClasses
  ) |>
    column_to_rownames("rownames")
  
  ## recode rowData factors with correct levels
  for(col in rowData_colSpecs$colName[rowData_colSpecs$colClass == "factor"]){
    rowData[[col]] <- factor(rowData[[col]], levels = strsplit(rowData_colSpecs$fctLevels[rowData_colSpecs$colName == col], ";|\\,\\ ")[[1]])
  }
  
  # metadata as list (if any)
  metadata.list <- tryCatch(
    jsonlite::read_json(file.path(tse.dir, "metadata.json"), simplifyVector = TRUE), 
    error = function(e) return(NULL), 
    warning = function(w) return(NULL))
  
  # assay(s) file names
  files_assays <- list.files(file.path(tse.dir, "assays"), pattern = "*.tsv", full.names = TRUE)
  
  assays <- lapply(files_assays,
                   fread) |>
    lapply(column_to_rownames, "rownames") |>
    lapply(as.matrix)
  names(assays) <- gsub(".tsv", "", basename(files_assays), fixed = TRUE)
  
  assaysOrder <- readLines(file.path(tse.dir, "assays", "assaysOrder.txt"))
  assays <- assays[assaysOrder]
  
  # read tree(s)
  rowTree <- tryCatch(
    tidytree::read.tree(file.path(tse.dir, "rowTree.tre")),
    error = function(e)
      return(NULL), warning = function(w)
        return(NULL)
  )
  colTree <- tryCatch(
    tidytree::read.tree(file.path(tse.dir, "colTree.tre")),
    error = function(e)
      return(NULL), warning = function(w)
        return(NULL)
  )
  
  # ensure correct row and columns order, exclude rows and columns that are not in common between
  # objects that should be in common
  # maybe this should actually shoot an error, because this is not supposed to 
  # happen
  assays_and_rowData <- append(lapply(assays, rownames), list("rowData" = rownames(rowData)))
  
  rows_in_common <- purrr::reduce(assays_and_rowData, intersect)
  
  assays_and_colData <- append(lapply(assays, colnames), list("rowData" = rownames(colData)))
  
  cols_in_common <- purrr::reduce(assays_and_colData, intersect)
  
  tse_rebuilt <- TreeSummarizedExperiment(
    metadata = metadata.list,
    assays = lapply(assays, function(x)
      x[rows_in_common, cols_in_common, drop = FALSE]),
    rowData = DataFrame(rowData[rows_in_common, , drop = FALSE]),
    colData = DataFrame(colData[cols_in_common, , drop = FALSE]),
    rowTree = rowTree,
    colTree = colTree
  )
  
  return(tse_rebuilt)
}


#' Read a TreeSummarizedExperiment from human readable files
#' 
#' This function is the reading counterpart of the `write_TSE_to_dir` function.
#' All files are human readable and maybe redundant, but at least more 
#' conservative.
#'
#' @param tse.dir \code{character}. The path that a `TreeSummarizedExperiment` 
#' was saved into
#'
#' @returns A \code{TreeSummarizedExperiment}
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#' 
#' # basic assay
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#' 
#' # sample data
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5), 
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#' 
#' # create rowData
#' rowData <- DataFrame(
#' Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#' Feature_gr2 = rep(c(2,3,4,6,1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#' 
#' # base tse
#' tse <- TreeSummarizedExperiment(
#'   assays = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#' 
#' # add altExp
#' alt_colData <- sampC
#' alt_colData$altTreat <- rep(0:1, 5)
#' 
#' alt_counts <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(assays = list(counts = alt_counts), rowData = alt_rowData, colData = alt_colData)
#' 
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#' 
#' tmpdir <- file.path(tempdir(), "testTSE_complete2")
#' unlink(tmpdir, recursive = TRUE)
#' write_TSE_to_dir(tse, tmpdir)
#' list.files(tmpdir, recursive = TRUE)
#' 
#' # re-load tse
#' tse2 <- read_TSE_from_dir(tmpdir)
#' all.equal(tse, tse2)  # TRUE

read_TSE_from_dir <- function(tse.dir) {
  tse_rebuilt <- read_TSE_from_dir_noAltExp(tse.dir)
  
  if(dir.exists(file.path(tse.dir, "altExps"))){
    # get names of the alternate expressions
    altExpNames <- readLines(file.path(tse.dir, "altExps", "AltExpOrder.txt"))
    
    for(x in altExpNames){
      tse_alt <- read_TSE_from_dir_noAltExp(file.path(tse.dir, "altExps", x))
      altExp(tse_rebuilt, x) <- tse_alt
    }
  }
  
  return(tse_rebuilt)
  
}