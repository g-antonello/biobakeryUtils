#' Specify classes of column data
#'
#' Captures the class of each column in a data frame or DataFrame, and for
#' factors, stores the levels as a JSON array string. This avoids separator
#' collision issues that arise when factor levels themselves contain delimiter
#' characters.
#'
#' @param df \code{data.frame} or \code{DataFrame} object
#'
#' @returns A \code{data.frame} with columns:
#'   \describe{
#'     \item{colName}{Column name, with \code{"rownames"} prepended}
#'     \item{colClass}{R class of the column}
#'     \item{fctLevels}{JSON array of factor levels, or \code{NA} for
#'       non-factors}
#'   }
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#'
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#'
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5),
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#'
#' rowData <- DataFrame(
#'   Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#'   Feature_gr2 = rep(c(2, 3, 4, 6, 1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#'
#' tse <- TreeSummarizedExperiment(
#'   assays  = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#'
#' DataFrameSpecs(colData(tse))
#' DataFrameSpecs(rowData(tse))

DataFrameSpecs <- function(df) {
  data.frame(
    colName = c("rownames", colnames(df)),
    colClass = c(
      "character",
      vapply(df, function(x) class(x)[1], character(1))
    ),
    # Store levels as a JSON array so any character can appear in a level
    # without risk of colliding with the delimiter.
    fctLevels = c(
      NA_character_,
      vapply(df, function(x) {
        if (is.factor(x)) {
          paste(levels(x), collapse = ";")
        } else {
          NA_character_
        }
      }, character(1))
    ),
    stringsAsFactors = FALSE
  )
}


#' Write a \code{TreeSummarizedExperiment} to a directory (no altExps)
#'
#' Serialises all slots of a \code{TreeSummarizedExperiment} to human-readable
#' files inside \code{out.dir}. The \code{altExp} slot is intentionally
#' ignored; use \code{\link{write_TSE_to_dir}} if you need altExps saved.
#'
#' The function errors if \code{out.dir} already exists, so it never silently
#' overwrites existing data.
#'
#' @param tse   A \code{TreeSummarizedExperiment}.
#' @param out.dir \code{character}. Path of the directory to create and
#'   populate.
#'
#' @importFrom data.table fwrite
#' @importFrom jsonlite write_json
#' @importFrom ape write.tree
#' @importFrom S4Vectors metadata
#' @importFrom TreeSummarizedExperiment rowTree colTree
#' @importFrom SummarizedExperiment colData rowData assay assayNames
#'
#' @returns \code{NULL} invisibly. Called for its side-effect of writing files.
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#'
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#'
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5),
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#'
#' rowData <- DataFrame(
#'   Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#'   Feature_gr2 = rep(c(2, 3, 4, 6, 1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#'
#' tse <- TreeSummarizedExperiment(
#'   assays  = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#'
#' tmpdir <- file.path(tempdir(), "testTSE_noAltExp")
#' write_TSE_to_dir_noAltExp(tse, tmpdir)
#' list.files(tmpdir, recursive = TRUE)
#' unlink(tmpdir, recursive = TRUE)

write_TSE_to_dir_noAltExp <- function(tse, out.dir) {
  
  if (dir.exists(out.dir)) {
    stop(
      "Output directory already exists: ", out.dir, "\n",
      "Remove it first or choose a different path."
    )
  }
  
  dir.create(out.dir, recursive = TRUE)
  
  # ---- colData ----------------------------------------------------------------
  colData.df <- as.data.frame(colData(tse))
  colData.df <- cbind(data.frame("rownames" = rownames(colData.df)), colData.df)
  
  data.table::fwrite(colData.df, file = file.path(out.dir, "colData.tsv"), sep = "\t")
  data.table::fwrite(
    DataFrameSpecs(colData(tse)),
    file = file.path(out.dir, "colData_colSpecs.tsv"),
    sep = "\t"
  )
  
  # ---- rowData ----------------------------------------------------------------
  rowData.df <- as.data.frame(rowData(tse)) 
  rowData.df <- cbind(data.frame("rownames" = rownames(rowData.df)), rowData.df)
  
  data.table::fwrite(rowData.df, file = file.path(out.dir, "rowData.tsv"), sep = "\t")
  data.table::fwrite(
    DataFrameSpecs(rowData(tse)),
    file = file.path(out.dir, "rowData_colSpecs.tsv"),
    sep = "\t"
  )
  
  # ---- metadata ---------------------------------------------------------------
  jsonlite::write_json(metadata(tse), file.path(out.dir, "metadata.json"))
  
  # ---- assays -----------------------------------------------------------------
  assays.dir <- file.path(out.dir, "assays")
  dir.create(assays.dir, recursive = TRUE)
  
  for (assayX in assayNames(tse)) {
    assay.df <- as.data.frame(assay(tse, assayX))
    assay.df <- cbind(data.frame(rownames = rownames(assay.df)), assay.df)
    
    data.table::fwrite(
      assay.df,
      file = file.path(assays.dir, paste0(assayX, ".tsv")),
      sep = "\t"
    )
  }
  # Record insertion order so reading can restore it exactly
  writeLines(assayNames(tse), con = file.path(assays.dir, "assaysOrder.txt"))
  
  # ---- trees ------------------------------------------------------------------
  if (inherits(tse, "TreeSummarizedExperiment")) {
    if (!is.null(rowTree(tse))) {
      ape::write.tree(rowTree(tse), file = file.path(out.dir, "rowTree.tre"))
    }
    if (!is.null(colTree(tse))) {
      # Bug fix: was incorrectly writing rowTree(tse) for both files
      ape::write.tree(colTree(tse), file = file.path(out.dir, "colTree.tre"))
    }
  }
  
  invisible(NULL)
}


#' Write a \code{TreeSummarizedExperiment} to a directory (including altExps)
#'
#' Serialises all slots of a \code{TreeSummarizedExperiment}, including any
#' \code{altExp} entries, to human-readable files. Each \code{altExp} is
#' written into a subdirectory \code{altExps/<name>/} using the same format.
#'
#' @param tse     A \code{TreeSummarizedExperiment}.
#' @param out.dir \code{character}. Path of the directory to create.
#'   Errors if the path already exists.
#'
#' @importFrom SingleCellExperiment altExpNames altExp
#'
#' @returns \code{NULL} invisibly. Called for its side-effect of writing files.
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#'
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#'
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5),
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#'
#' rowData <- DataFrame(
#'   Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#'   Feature_gr2 = rep(c(2, 3, 4, 6, 1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#'
#' tse <- TreeSummarizedExperiment(
#'   assays  = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#'
#' alt_counts  <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(
#'   assays  = list(counts = alt_counts),
#'   rowData = alt_rowData,
#'   colData = sampC
#' )
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#'
#' tmpdir <- file.path(tempdir(), "testTSE_complete")
#' write_TSE_to_dir(tse, tmpdir)
#' list.files(tmpdir, recursive = TRUE)
#' unlink(tmpdir, recursive = TRUE)

write_TSE_to_dir <- function(tse, out.dir) {
  
  # out.dir existence is already checked inside write_TSE_to_dir_noAltExp,
  # but we check here first so the error message is unambiguous about the
  # top-level path before we attempt to create any subdirectories.
  if (dir.exists(out.dir)) {
    stop(
      "Output directory already exists: ", out.dir, "\n",
      "Remove it first or choose a different path."
    )
  }
  
  write_TSE_to_dir_noAltExp(tse = tse, out.dir = out.dir)
  
  if (length(altExpNames(tse)) > 0) {
    for (aexp in altExpNames(tse)) {
      altexp.outdir <- file.path(out.dir, "altExps", aexp)
      # write_TSE_to_dir_noAltExp will create and guard altexp.outdir itself
      write_TSE_to_dir_noAltExp(altExp(tse, aexp), out.dir = altexp.outdir)
    }
    writeLines(
      altExpNames(tse),
      con = file.path(out.dir, "altExps", "AltExpOrder.txt")
    )
  }
  
  invisible(NULL)
}


#' Read a \code{TreeSummarizedExperiment} from a directory (no altExps)
#'
#' Counterpart of \code{\link{write_TSE_to_dir_noAltExp}}. Reconstructs all
#' slots except \code{altExp}. Factor levels are restored exactly, including
#' unused levels, because they were stored as JSON arrays.
#'
#' @param tse.dir \code{character}. Path written by
#'   \code{write_TSE_to_dir_noAltExp} or \code{write_TSE_to_dir}.
#'
#' @importFrom data.table fread set
#' @importFrom ape read.tree
#' @importFrom jsonlite read_json
#' @importFrom S4Vectors DataFrame
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @returns A \code{TreeSummarizedExperiment}.
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#'
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#'
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5),
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#'
#' rowData <- DataFrame(
#'   Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#'   Feature_gr2 = rep(c(2, 3, 4, 6, 1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#'
#' tse <- TreeSummarizedExperiment(
#'   assays  = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#' metadata(tse)$info <- "test_metadata"
#'
#' tmpdir <- file.path(tempdir(), "testTSE_read")
#' write_TSE_to_dir_noAltExp(tse, tmpdir)
#' tse2 <- read_TSE_from_dir_noAltExp(tmpdir)
#' tse2
#' unlink(tmpdir, recursive = TRUE)

read_TSE_from_dir_noAltExp <- function(tse.dir) {
  
  # ---- helper: apply colSpecs to a freshly-read data.table ------------------
  .apply_specs <- function(dt, specs) {
    # dt has already been read with colClasses applied; we only need to
    # re-apply factor levels (which fread cannot do).
    factor_rows <- which(specs$colClass == "factor")
    for (i in factor_rows) {
      lvls <- strsplit(unlist(specs[i, "fctLevels"]), ";")[[1]]
      data.table::set(dt, j = i, value = factor(dt[[i]], levels = lvls))
    }
    return(dt)
  }
  
  # ---- colData ---------------------------------------------------------------
  colData_specs <- data.table::fread(
    file.path(tse.dir, "colData_colSpecs.tsv"), sep = "\t", header = TRUE
  )
  colClasses_col        <- colData_specs$colClass
  names(colClasses_col) <- colData_specs$colName
  
  colData.df <- data.table::fread(
    file.path(tse.dir, "colData.tsv"), sep = "\t", colClasses = colClasses_col, 
    data.table = FALSE
  )
  colData.df <- .apply_specs(dt = colData.df, specs = colData_specs)
  rownames(colData.df) <- colData.df$rownames
  colData.df$rownames <- NULL
  
  # ---- rowData ---------------------------------------------------------------
  rowData_specs <- data.table::fread(
    file.path(tse.dir, "rowData_colSpecs.tsv"), sep = "\t", header = TRUE
  )
  colClasses_row        <- rowData_specs$colClass
  names(colClasses_row) <- rowData_specs$colName
  
  rowData.df <- data.table::fread(
    file.path(tse.dir, "rowData.tsv"), sep = "\t", colClasses = colClasses_row,
    data.table = FALSE
  )
  rowData.df <- .apply_specs(rowData.df, rowData_specs)
  rownames(rowData.df) <- rowData.df$rownames
  rowData.df$rownames <- NULL
  
  # ---- metadata --------------------------------------------------------------
  metadata.path <- file.path(tse.dir, "metadata.json")
  metadata.list <- if (!file.exists(metadata.path)) {
    NULL  # no metadata saved — perfectly normal, no message needed
  } else {
    tryCatch(
      read_json(metadata.path, simplifyVector = TRUE),
      error   = function(e) {
        message("metadata.json exists but could not be parsed: ", e$message,
                "\nProceeding with empty metadata.")
        NULL
      },
      warning = function(w) {
        message("Warning parsing metadata.json: ", w$message,
                "\nProceeding with empty metadata.")
        NULL
      }
    )
  }
  
  # ---- assays ----------------------------------------------------------------
  assays.dir   <- file.path(tse.dir, "assays")
  assaysOrder  <- readLines(file.path(assays.dir, "assaysOrder.txt"))
  
  # Use an explicit regex so assaysOrder.txt is never accidentally picked up
  files_assays <- list.files(assays.dir, pattern = "\\.tsv$", full.names = TRUE)
  assays <- lapply(files_assays, function(f) {
    tmp <- data.table::fread(f, data.table = FALSE)
    rownames(tmp) <- tmp$rownames
    tmp$rownames <- NULL
      return(as.matrix(tmp))
  })
  names(assays) <- gsub("\\.tsv$", "", basename(files_assays))
  
  # Restore original insertion order
  assays <- assays[assaysOrder]
  
  # ---- trees -----------------------------------------------------------------
  .read_tree_safe <- function(path) {
    if (!file.exists(path)) return(NULL)  # absent tree is normal — no message
    tryCatch(
      ape::read.tree(path),
      error   = function(e) {
        message("Tree file '", path, "' exists but could not be read: ", e$message,
                "\nProceeding without this tree.")
        NULL
      },
      warning = function(w) {
        message("Warning reading tree file '", path, "': ", w$message,
                "\nProceeding without this tree.")
        NULL
      }
    )
  }
  rowTree <- .read_tree_safe(file.path(tse.dir, "rowTree.tre"))
  colTree <- .read_tree_safe(file.path(tse.dir, "colTree.tre"))
  
  # ---- dimension consistency check -------------------------------------------
  # Mismatches should never occur for files written by this package, but if
  # they do we error loudly rather than silently dropping rows/columns.
  assay_rownames <- lapply(assays, rownames)
  assay_colnames <- lapply(assays, colnames)
  
  for (nm in names(assay_rownames)) {
    if (!identical(sort(assay_rownames[[nm]]), sort(rownames(rowData.df)))) {
      stop(
        "Row names of assay '", nm, "' do not match rowData row names.\n",
        "This indicates a corrupted or manually edited directory."
      )
    }
    if (!identical(sort(assay_colnames[[nm]]), sort(rownames(colData.df)))) {
      stop(
        "Column names of assay '", nm, "' do not match colData row names.\n",
        "This indicates a corrupted or manually edited directory."
      )
    }
  }
  
  # Reindex assays to rowData / colData order (which is the authoritative order)
  row_order <- rownames(rowData.df)
  col_order <- rownames(colData.df)
  assays <- lapply(assays, function(x) x[row_order, col_order, drop = FALSE])
  
  # ---- reconstruct TSE -------------------------------------------------------
  TreeSummarizedExperiment(
    metadata = metadata.list,
    assays   = assays,
    rowData  = S4Vectors::DataFrame(rowData.df),
    colData  = S4Vectors::DataFrame(colData.df),
    rowTree  = rowTree,
    colTree  = colTree
  )
}


#' Read a \code{TreeSummarizedExperiment} from a directory (including altExps)
#'
#' Counterpart of \code{\link{write_TSE_to_dir}}. Reads the main experiment
#' and all \code{altExp} entries back from disk.
#'
#' @param tse.dir \code{character}. Path written by \code{write_TSE_to_dir}.
#'
#' @importFrom SingleCellExperiment altExp<-
#'
#' @returns A \code{TreeSummarizedExperiment}.
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' set.seed(1234)
#'
#' count.mat <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count.mat) <- tinyTree$tip.label
#' colnames(count.mat) <- paste0("C_", 1:10)
#'
#' sampC <- data.frame(
#'   condition = factor(rep(c("control", "trt"), each = 5),
#'                      levels = c("trt", "control")),
#'   gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
#'                   levels = c("male", "female"))
#' )
#' rownames(sampC) <- colnames(count.mat)
#'
#' rowData <- DataFrame(
#'   Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
#'   Feature_gr2 = rep(c(2, 3, 4, 6, 1), 2)
#' )
#' rownames(rowData) <- rownames(count.mat)
#'
#' tse <- TreeSummarizedExperiment(
#'   assays  = list(counts = count.mat),
#'   colData = sampC,
#'   rowData = rowData,
#'   rowTree = tinyTree
#' )
#'
#' alt_counts  <- count.mat[1:5, ]
#' alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
#' alt_se <- TreeSummarizedExperiment(
#'   assays  = list(counts = alt_counts),
#'   rowData = alt_rowData,
#'   colData = sampC
#' )
#' altExp(tse, "altSubset") <- alt_se
#' metadata(tse)$info <- "test_metadata"
#'
#' tmpdir <- file.path(tempdir(), "testTSE_complete2")
#' write_TSE_to_dir(tse, tmpdir)
#' tse2 <- read_TSE_from_dir(tmpdir)
#' all.equal(tse, tse2)
#' unlink(tmpdir, recursive = TRUE)

read_TSE_from_dir <- function(tse.dir) {
  tse_rebuilt <- read_TSE_from_dir_noAltExp(tse.dir)
  
  altExps.dir <- file.path(tse.dir, "altExps")
  if (dir.exists(altExps.dir)) {
    altExp_names <- readLines(file.path(altExps.dir, "AltExpOrder.txt"))
    for (x in altExp_names) {
      tse_alt <- read_TSE_from_dir_noAltExp(file.path(altExps.dir, x))
      altExp(tse_rebuilt, x) <- tse_alt
    }
  }
  
  tse_rebuilt
}