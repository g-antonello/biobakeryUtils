#' Specify classes of column data
#' Most importantly, the levels of a factor
#' 
#' @param tse \code{TreeSummarizedExperiment} Object
#'
#' @returns A \code{data.frame}
#' @export
#'
#' @examples
#'  library(TreeSummarizedExperiment)
#'  data("tinyTree")
#'
#'  # the count.mat table
#'  count.mat <- matrix(rpois(100, 50), nrow = 10)
#'
#'  rownames(count.mat) <- c(tinyTree$tip.label)
#'
#'  colnames(count.mat) <- paste("C_", 1:10, sep = "_")
#'
#'  # The sample information
#' sampC <- data.frame(
#'   condition = rep(c("control", "trt"), each = 5),
#'   gender = sample(x = 1:2, size = 10, replace = TRUE)
#'   )
#'
#' rownames(sampC) <- colnames(count.mat)
#'
#'  # build a TreeSummarizedExperiment object
#'  tse <- TreeSummarizedExperiment(
#'    assays = list("counts" = count.mat),
#'    colData = sampC,
#'    rowTree = tinyTree
#'    )
#'    
#'  tse

colDataSpecs <- function(tse){
  
  colSpecs.df <- data.frame(
    colName = colnames(colData(tse)),
    colClass = sapply(colData(tse), class),
    fctLevels = sapply(colData(tse), function(x) paste(levels(x), collapse = ", "))
  )
  
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
#' @importFrom readr write_tsv
#' @importFrom jsonlite write_json
#' @importFrom tibble rownames_to_column
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom ape write.tree
#' @importFrom S4Vectors metadata
#'
#' @returns Nothing, the main output of this function is multiple .tsv files
#'
#' @examples
#'  library(TreeSummarizedExperiment)
#'  data("tinyTree")
#'
#'  # the count.mat table
#'  count.mat <- matrix(rpois(100, 50), nrow = 10)
#'
#'  rownames(count.mat) <- c(tinyTree$tip.label)
#'
#'  colnames(count.mat) <- paste("C_", 1:10, sep = "_")
#'
#'  # The sample information
#' sampC <- data.frame(
#'   condition = rep(c("control", "trt"), each = 5),
#'   gender = sample(x = 1:2, size = 10, replace = TRUE)
#'   )
#'
#' rownames(sampC) <- colnames(count.mat)
#'
#'  # build a TreeSummarizedExperiment object
#'  tse <- TreeSummarizedExperiment(
#'    assays = list("counts" = count.mat),
#'    colData = sampC,
#'    rowTree = tinyTree
#'    )
#'
#'  # write the tse
#'
#'  set.seed(1234)
#'
#'  write_TSE_to_dir_noAltExp(tse, tempdir())
#'
#'  list.files(file.path(tempdir(), "tse"), recursive = TRUE)

write_TSE_to_dir_noAltExp <- function(tse, out.dir) {
  # write the colData
  write_tsv(
    colData(tse) |> as.data.frame() |> rownames_to_column("rownames"),
    file = file.path(out.dir, paste0("colData", ".tsv"))
  )
  # also write column specifications
  write_tsv(
    colDataSpecs(tse), file.path(out.dir, paste0("colData_colSpecs", ".tsv"))
  )
  
  # write the rowData
  write_tsv(
    rowData(tse) |> as.data.frame() |> rownames_to_column("rownames"),
    file = file.path(out.dir, paste0("rowData", ".tsv"))
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

write_TSE_to_dir <- function(tse, out.dir){
  if(dir.exists(out.dir)){
    error("Output directory already exists")
  } 
  dir.create(out.dir, recursive = TRUE)
  
  write_TSE_to_dir_noAltExp(tse = tse, out.dir = out.dir)
  
  if(length(altExpNames(tse)) > 0){
  for(aexp in altExpNames(tse)){
    
    altexp.outdir <- file.path(out.dir, "altExps", aexp)
    
    if(dir.exists(altexp.outdir)){
      error("Output directory already exists")
    } 
    dir.create(altexp.outdir, recursive = TRUE)
    
    write_TSE_to_dir_noAltExp(altExp(tse, aexp), out.dir = altexp.outdir)
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
#' @importFrom readr read_tsv
#' @importFrom tibble column_to_rownames
#' @importFrom purrr reduce
#' @importFrom jsonlite read_json
#'
#' @examples
#'  library(TreeSummarizedExperiment)
#'  data("tinyTree")
#'
#'  # the count.mat table
#'  count.mat <- matrix(rpois(100, 50), nrow = 10)
#'
#'  rownames(count.mat) <- c(tinyTree$tip.label)
#'
#'  colnames(count.mat) <- paste("C_", 1:10, sep = "_")
#'
#'  # The sample information
#' sampC <- data.frame(
#'   condition = rep(c("control", "trt"), each = 5),
#'   gender = as.numeric(sample(x = 1:2, size = 10, replace = TRUE))
#'   )
#'
#' rownames(sampC) <- colnames(count.mat)
#'
#'  # build a TreeSummarizedExperiment object
#'  tse <- TreeSummarizedExperiment(
#'    assays = list("counts" = count.mat),
#'    colData = sampC,
#'    rowTree = tinyTree
#'    )
#'  tse <- reorder_taxa_with_phyloTree_labels(tse)
#'  # save the tse
#'
#'  set.seed(1234)
#'
#'  write_TSE_to_dir(tse, tempdir())
#'
#'  list.files(file.path(tempdir(), "tse"), recursive = TRUE)
#'
#'  tse2 <- read_TSE_from_dir(file.path(tempdir(), "tse"))
#'
#'  identical(tse, tse2) # FALSE
#'  all.equal(tse, tse2) # TRUE
#'  identical(rownames(tse), rownames(tse2)) # TRUE
#'  identical(colnames(tse), colnames(tse2)) # TRUE
#'  identical(tse@rowTree$phylo$tip.label, tse2@rowTree$phylo$tip.label) # TRUE
#'  identical(tse@rowTree$phylo$node.label, tse2@rowTree$phylo$node.label) # TRUE
#'  identical(
#'  round(tse@rowTree$phylo$edge.length,10),
#'  round(tse2@rowTree$phylo$edge.length, 10)
#'  ) # TRUE, identical at least to the 10th digit
#'
#'  identical(assay(tse), assay(tse)) # FALSE, I don't yet know why
#'  identical(round(assay(tse), 10), round(assay(tse), 10)) # The encoding
#'  # could be improved, but the matrices are practically the same

read_TSE_from_dir_noAltExp <- function(tse.dir) {
  
  # colData
  ## get colClasses and specify them while reading colData
  colDataSpecs <- read_tsv(
    file.path(tse.dir, "colData_colSpecs.tsv"), show_col_types = FALSE) 
    # this is done here because read_tsv may change in the future
    colDataSpecs$readrClass <- ifelse(colDataSpecs$colClass == "factor", 
                                       "c", 
                                       substr(colDataSpecs$colClass, 1, 1)
                                       )
  
  ## read colData and specify colClasses
  colData <- read_tsv(
    file.path(tse.dir, "colData.tsv"),
    col_types = colDataSpecs$readrClass,
    progress = FALSE
  ) |>
    column_to_rownames("rownames")
  
  ## recode colData factors with correct levels
  for(col in colDataSpecs$colName[colDataSpecs$colClass == "factor"]){
    colData[[col]] <- factor(colData[[col]], levels = strsplit(colDataSpecs$fctLevels[colDataSpecs$colName == col], "\\,\\ ")[[1]])
  }
  
  # rowData
  rowData <- read_tsv(
    file.path(tse.dir, "rowData.tsv"),
    progress = FALSE,
    show_col_types = FALSE
  ) |>
    column_to_rownames("rownames")
  
  # metadata as list (if any)
  metadata.list <- tryCatch(
    jsonlite::read_json(file.path(tse.dir, "metadata.json"), simplifyVector = TRUE), 
    error = function(e) return(NULL), 
    warning = function(w) return(NULL))
  
  # assay(s) file names
  files_assays <- list.files(file.path(tse.dir, "assays"), pattern = "*.tsv", full.names = TRUE)
  
  assays <- lapply(files_assays,
                   read_tsv,
                   progress = FALSE,
                   show_col_types = FALSE) |>
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
    rowData = DataFrame(rowData[rows_in_common, ]),
    colData = DataFrame(colData[cols_in_common, ]),
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
#'  library(TreeSummarizedExperiment)
#'  data("tinyTree")
#'
#'  # the count.mat table
#'  count.mat <- matrix(rpois(100, 50), nrow = 10)
#'
#'  rownames(count.mat) <- c(tinyTree$tip.label)
#'
#'  colnames(count.mat) <- paste("C_", 1:10, sep = "_")
#'
#'  # The sample information
#' sampC <- data.frame(
#'   condition = rep(c("control", "trt"), each = 5),
#'   gender = as.numeric(sample(x = 1:2, size = 10, replace = TRUE))
#'   )
#'
#' rownames(sampC) <- colnames(count.mat)
#'
#'  # build a TreeSummarizedExperiment object
#'  tse <- TreeSummarizedExperiment(
#'    assays = list("counts" = count.mat),
#'    colData = sampC,
#'    rowTree = tinyTree
#'    )
#'  tse <- reorder_taxa_with_phyloTree_labels(tse)
#'  # save the tse
#'
#'  set.seed(1234)
#'
#'  write_TSE_to_dir(tse, tempdir())
#'
#'  list.files(file.path(tempdir(), "tse"), recursive = TRUE)
#'
#'  tse2 <- read_TSE_from_dir(file.path(tempdir(), "tse"))
#'
#'  identical(tse, tse2) # FALSE
#'  all.equal(tse, tse2) # TRUE
#'  identical(rownames(tse), rownames(tse2)) # TRUE
#'  identical(colnames(tse), colnames(tse2)) # TRUE
#'  identical(tse@rowTree$phylo$tip.label, tse2@rowTree$phylo$tip.label) # TRUE
#'  identical(tse@rowTree$phylo$node.label, tse2@rowTree$phylo$node.label) # TRUE
#'  identical(
#'  round(tse@rowTree$phylo$edge.length,10),
#'  round(tse2@rowTree$phylo$edge.length, 10)
#'  ) # TRUE, identical at least to the 10th digit
#'
#'  identical(assay(tse), assay(tse)) # FALSE, I don't yet know why
#'  identical(round(assay(tse), 10), round(assay(tse), 10)) # The encoding
#'  # could be improved, but the matrices are practically the same

read_TSE_from_dir <- function(tse.dir) {
  tse_rebuilt <- read_TSE_from_dir_noAltExp(tse.dir)
  
  if(dir.exists(file.path(tse.dir, "altExps"))){
    altexps_dirs <- list.dirs(file.path(tse.dir, "altExps"), recursive = FALSE)
    
    altExprs <- lapply(altexps_dirs, function(altexp.dir) read_TSE_from_dir_noAltExp(altexp.dir))
    
    exps_order <- readLines(file.path(tse.dir, "altExps", "AltExpOrder.txt"))
    
    altExprs <- altExprs[exps_order]
    altExp(tse_rebuilt) <- altExprs
  }
  
  
  return(tse_rebuilt)
  }