

#' Split a `TreeSummarizedExperiment` object into multiple .tsv files
#'
#' The components are saved into a directory that gets created by this function
#'
#' @param tse \code{TreeSummarizedExperiment} object with or without
#' rowTree and/or colTree
#' @param dir \code{character} specifying the directory in which a new directory
#' will be generated, with the name of the R object passed to the `tse`
#' parameter.
#'
#' @importFrom readr write_tsv
#' @importFrom tibble rownames_to_column
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom ape write.tree
#'
#' @returns Nothing, the main output of this function is multiple .tsv files
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
#'  # write the tse
#'
#'  set.seed(1234)
#'
#'  write_TSE_to_dir(tse, tempdir())
#'
#'  list.files(file.path(tempdir(), "tse"), recursive = TRUE)

write_TSE_to_dir <- function(tse, dir) {
  # create output directory if it doesn't exist
  tse_dirName <- deparse(substitute(tse))
  dir.create(file.path(dir, tse_dirName),
             showWarnings = FALSE,
             recursive = TRUE)
  
  # write the colData
  write_tsv(
    colData(tse) |> as.data.frame() |> rownames_to_column("rownames"),
    file = file.path(dir, tse_dirName, paste0("colData", ".tsv"))
  )
  
  # write the rowData
  write_tsv(
    rowData(tse) |> as.data.frame() |> rownames_to_column("rownames"),
    file = file.path(dir, tse_dirName, paste0("rowData", ".tsv"))
  )
  
  # write the assays, including the subdirectory
  dir.create(
    file.path(dir, tse_dirName, "assays"),
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  for (assayX in assayNames(tse)) {
    write_tsv(
      assay(tse, assayX) |>
        as.data.frame() |>
        rownames_to_column("rownames"),
      file = file.path(dir, tse_dirName, "assays", paste0(assayX, ".tsv"))
    )
    
  }
  
  # write rowTree file
  if (inherits(tse, "TreeSummarizedExperiment")) {
    if (!is.null(rowTree(tse))) {
      write.tree(rowTree(tse), file = file.path(dir, tse_dirName, "rowTree.tre"))
    }
    # write colTree file
    if (!is.null(colTree(tse))) {
      write.tree(rowTree(tse), file = file.path(dir, tse_dirName, "colTree.tre"))
    }
  }
}

#' read multiple tsv files into a TreeSummarizedExperiment
#'
#' This utility function is the counterpart of
#'
#' @param dir \code{character} indicating the directory that was created
#' with `write_TSE_to_dir()`
#'
#' @returns A \code{TreeSummarizedExperiment}
#'
#' @importFrom readr read_tsv
#' @importFrom tibble column_to_rownames
#' @importFrom purrr reduce
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
#'  tse2 <- read_TSE_from_dir(dir = file.path(tempdir(), "tse"))
#'
#'  identical(tse, tse2) # FALSE
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

read_TSE_from_dir <- function(dir) {
  files <- list.files(dir, full.names = TRUE)
  colData <- read_tsv(
    grep("colData", files, value = TRUE),
    progress = FALSE,
    show_col_types = FALSE
  ) |>
    column_to_rownames("rownames")
  rowData <- read_tsv(
    grep("rowData", files, value = TRUE),
    progress = FALSE,
    show_col_types = FALSE
  ) |>
    column_to_rownames("rownames")
  
  # read assay file names
  files_assays <- list.files(file.path(dir, "assays"), full.names = TRUE)
  
  assays <- lapply(files_assays,
                   read_tsv,
                   progress = FALSE,
                   show_col_types = FALSE) |>
    lapply(column_to_rownames, "rownames") |>
    lapply(as.matrix)
  names(assays) <- gsub(".tsv", "", basename(files_assays), fixed = TRUE)
  
  # read tree(s)
  rowTree <- tryCatch(
    tidytree::read.tree(grep("rowTree.tre", files, value = TRUE)),
    error = function(e)
      return(NULL)
  )
  colTree <- tryCatch(
    tidytree::read.tree(grep("colTree.tre", files, value = TRUE)),
    error = function(e)
      return(NULL)
  )
  
  # ensure correct row and columns order, exclude rows and columns that are not in common between
  # objects that should be in common
  assays_and_rowData <- append(lapply(assays, rownames), list("rowData" = rownames(rowData)))
  
  rows_in_common <- purrr::reduce(assays_and_rowData, intersect)
  
  assays_and_colData <- append(lapply(assays, colnames), list("rowData" = rownames(colData)))
  
  cols_in_common <- purrr::reduce(assays_and_colData, intersect)
  
  tse_rebuilt <- TreeSummarizedExperiment(
    assays = lapply(assays, function(x)
      x[rows_in_common, cols_in_common]),
    rowData = DataFrame(rowData[rows_in_common, ]),
    colData = DataFrame(colData[cols_in_common, ]),
    rowTree = rowTree,
    colTree = colTree
  )
  
  # No need, but I will leave it here to be sure
  tse_rebuilt_reordered_rows <- reorder_taxa_with_phyloTree_labels(tse_rebuilt)
  
  return(tse_rebuilt_reordered_rows)
}
