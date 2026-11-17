make_test_tse <- function() {
  set.seed(1234)
  suppressMessages(library(TreeSummarizedExperiment))
  suppressMessages(library(treeio))
  data("tinyTree", package = "TreeSummarizedExperiment")
  #-------------------
  # force tinyTree to have canonical ape/treeio elements order.
  # When non-example trees are read, no errors are given because newick tree
  # are read and stored in a canonical way. this one would fail some tests
  tinyTree_backup <- tinyTree
  canonical <- c("edge", "edge.length", "Nnode", "node.label", "tip.label")
  presentNames <- canonical[canonical %in% names(tinyTree)]
  restNames <- setdiff(names(tinyTree), presentNames)
  tinyTree <- tinyTree[c(presentNames, restNames)]  # reorder in-place, preserves class
  tinyTree <- structure(unclass(tinyTree)[c(presentNames, restNames)], class = "phylo")
  attr(tinyTree, "order") <- "cladewise"
  #-------------------
  if (!all.equal(tinyTree_backup, tinyTree)){
    stop("tree conversion not succesful while making example TSE object")
  }
  
  rowTree <- tinyTree
  colTree <- tinyTree
  colTree$tip.label <- gsub("t", "C_", colTree$tip.label)
  
  # basic assay
  count.mat <- matrix(rpois(100, 50), nrow = 10)
  rownames(count.mat) <- tinyTree$tip.label
  colnames(count.mat) <- paste0("C_", 1:10)
  
  # colData
  sampC <- data.frame(
    condition = factor(rep(c("control", "trt"), each = 5), levels = c("trt", "control")),
    gender = factor(sample(
      c("male", "female"), size = 10, replace = TRUE
    ), levels = c("male", "female")),
    case_withNAs = c(1.21, 0.87, 1.22, NA, 0.03, -0.11, -1.87, NA, 15, NA)
  )
  rownames(sampC) <- colnames(count.mat)
  
  # create rowData
  rowData <- DataFrame(Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
                       Feature_gr2 = rep(c(2, 3, 4, 6, 1), 2))
  rownames(rowData) <- rownames(count.mat)
  
  # base tse
  tse <- TreeSummarizedExperiment(
    assays = list(counts = count.mat),
    colData = sampC,
    rowData = rowData,
    rowTree = rowTree,
    colTree = colTree
  )
  
  # add altExp
  alt_colData <- sampC
  alt_colData$altTreat <- rep(0:1, 5)
  
  alt_counts <- count.mat[1:5, ]
  alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
  alt_se <- TreeSummarizedExperiment(
    assays = list(counts = alt_counts),
    rowData = alt_rowData,
    colData = alt_colData
  )
  
  altExp(tse, "altSubset") <- alt_se
  
  # add metadata
  metadata(tse) <- list("test_chr" = "test_metadata", "test_df" = mtcars)
  
  return(tse)
}

# --------------------- DataFrameSpecs tests ----------------------------------

test_that("DataFrameSpecs returns correct structure", {
  
  # tests on coldata
  tse <- make_test_tse()
  specs <- DataFrameSpecs(colData(tse))
  
  expect_s3_class(specs, "data.frame")
  expect_true(all(c("colName", "colClass", "fctLevels") %in% colnames(specs)))
  
  # condition should be factor with correct levels
  cond_row <- specs[specs$colName == "condition", ]
  expect_equal(cond_row$colClass, "factor")
  expect_true(all(c("trt", "control") %in% strsplit(cond_row$fctLevels, ";")[[1]]))
  
  # same test on rowData
  specs <- DataFrameSpecs(rowData(tse))
  
  expect_s3_class(specs, "data.frame")
  expect_true(all(c("colName", "colClass", "fctLevels") %in% colnames(specs)))
  
  # Feature_gr1 should be factor with correct levels
  cond_row <- specs[specs$colName == "Feature_gr1", ]
  expect_equal(cond_row$colClass, "factor")
  expect_true(all(c("blue", "red") %in% strsplit(cond_row$fctLevels, ";")[[1]]))
})
# -----------------------------------------------------------------------------

# -------  test I/O of the _noAltExp() version of function  -------------------

test_that("TSE round-trips without altExp", {
  tse <- make_test_tse()
  
  tmpDestDir <- "testOutDir"
  unlink(tmpDestDir, recursive = TRUE)
  write_TSE_to_dir_noAltExp(tse, out.dir = tmpDestDir)
  expect_error(write_TSE_to_dir_noAltExp(tse, out.dir = tmpDestDir))  # shoot error if the target directory exists
  
  # expect existence of the following data
  expect_true(file.exists(file.path(tmpDestDir, "colData.tsv")))
  expect_true(file.exists(file.path(tmpDestDir, "rowData.tsv")))
  expect_true(file.exists(file.path(tmpDestDir, "rowTree.tre")))
  expect_true(file.exists(file.path(tmpDestDir, "colTree.tre")))
  
  # expect that an altExp does not exist, because this is the _noAltExp type
  # of funciton
  expect_false(dir.exists(file.path(tmpDestDir, "altExps")))
  
  tse2 <- read_TSE_from_dir_noAltExp(tse.dir = tmpDestDir)
  
  # metadata, assays, row/col order preserved
  expect_equal(metadata(tse)$info, metadata(tse2)$info)
  expect_equal(rownames(tse), rownames(tse2))
  expect_equal(colnames(tse), colnames(tse2))
  expect_equal(dim(assay(tse)), dim(assay(tse2)))
  
  # factors preserved
  expect_true(is.factor(colData(tse2)$condition))
  expect_equal(levels(colData(tse)$condition),
               levels(colData(tse2)$condition))
  
  # expect rowLinks to be different, but column by column being identical
  expect_equal(rowLinks(tse), rowLinks(tse2))
  expect_equal(rowLinks(tse), rowLinks(tse2))
  expect_true(all.equal(rowLinks(tse), rowLinks(tse2)))
  
  # expect that AltExps are not saved
  expect_true(length(altExpNames(tse)) > 0 & length(altExpNames(tse2)) == 0)
  
  # last, and most important check, check if they are ultimately equal overall 
  # but ignoring the altExps
  altExp(tse) <- NULL
  expect_equal(tse, tse2)
  
  # clean up before next test
  unlink(tmpDestDir, recursive = TRUE)
})

# -----------------------------------------------------------------------------

# test full I/O functions (including AltExp)

test_that("TSE round-trips with altExp", {
  tse <- make_test_tse()
  tmpDest <- file.path(tempdir(), "testTSE")
  write_TSE_to_dir(tse, tmpDest)
  
  # expect existence of a bunch of essential files for building a TSE
  expect_true(file.exists(file.path(tmpDest, "colData.tsv")))
  expect_true(file.exists(file.path(tmpDest, "rowData.tsv")))
  expect_true(file.exists(file.path(tmpDest, "rowTree.tre")))
  expect_true(file.exists(file.path(tmpDest, "colTree.tre")))
  # in contrast with the example above, now we expect a TRUE here below
  expect_true(dir.exists(file.path(tmpDest, "altExps")))
  
  tse2 <- read_TSE_from_dir(tmpDest)
  
  # check altExp presence
  expect_true("altSubset" %in% altExpNames(tse2))
  expect_equal(rownames(altExp(tse, "altSubset")), 
               rownames(altExp(tse2, "altSubset")))
  
  # check metadata persisted
  expect_equal(metadata(tse)$info, metadata(tse2)$info)
  
  # last, and most important check, check if they are ultimately equal overall
  expect_equal(tse, tse2)
  
  # clean up before next test
  unlink(tmpDest, recursive = TRUE)
})

# Expected errors in blank directory rules

test_that("write_TSE_to_dir detects existing directory", {
  tse <- make_test_tse()
  
  tmpDest <- "testTSE"
  
  # expect no error on the _noAltExp version of the code for the first attempt
  expect_no_error(write_TSE_to_dir_noAltExp(tse, tmpDest))
  # error if writing is re-attempted 
  expect_error(write_TSE_to_dir_noAltExp(tse, tmpDest))
  
  unlink(tmpDest, recursive = TRUE)
  
  # expect no error on the full code
  expect_no_error(write_TSE_to_dir(tse, tmpDest))
  # expect error if writing is re-attempted
  expect_error(write_TSE_to_dir(tse, tmpDest))
  
  # expect an error if the directory already exists, 
  # even if it's empty
  tmpDest2 <- file.path("testTSE2")
  dir.create(tmpDest2)
  expect_error(write_TSE_to_dir(tse, tmpDest2))
  
  # clean up before next test
  unlink(tmpDest, recursive = TRUE)
  unlink(tmpDest2, recursive = TRUE)
})


# test whether extra levels of the factor get preserved even if they are not 
# in the colData anymore

test_that("Factors with unused levels are preserved", {
  tse <- make_test_tse()
  tse$condition <- factor(tse$condition, levels = c("trt", "control", "extra"))
  
  tmpdir <- tempfile("testTSE_factors")
  write_TSE_to_dir(tse, tmpdir)
  tse2 <- read_TSE_from_dir_noAltExp(tmpdir)
  
  expect_equal(levels(tse$condition), levels(tse2$condition))
})

test_that("colTree round-trips correctly", {
  tse <- make_test_tse()
  tmpDest <- file.path(tempdir(), "testTSE_colTree")
  write_TSE_to_dir(tse, tmpDest)
  tse2 <- read_TSE_from_dir(tmpDest)
  expect_equal(colTree(tse)$tip.label, colTree(tse2)$tip.label)
  unlink(tmpDest, recursive = TRUE)
})

test_that("multiple assays round-trip in correct order", {
  tse <- make_test_tse()
  assay(tse, "logcounts") <- log1p(assay(tse, "counts"))
  tmpDest <- file.path(tempdir(), "testTSE_assays")
  write_TSE_to_dir(tse, tmpDest)
  tse2 <- read_TSE_from_dir(tmpDest)
  expect_equal(assayNames(tse), assayNames(tse2))
  unlink(tmpDest, recursive = TRUE)
})
