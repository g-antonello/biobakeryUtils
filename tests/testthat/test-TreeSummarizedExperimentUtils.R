make_test_tse <- function() {
library(TreeSummarizedExperiment)
data("tinyTree", package = "TreeSummarizedExperiment")
set.seed(1234)

# basic assay
count.mat <- matrix(rpois(100, 50), nrow = 10)
rownames(count.mat) <- tinyTree$tip.label
colnames(count.mat) <- paste0("C_", 1:10)

# sample data
sampC <- data.frame(
  condition = factor(rep(c("control", "trt"), each = 5), 
                     levels = c("trt", "control")),
  gender = factor(sample(c("male", "female"), size = 10, replace = TRUE),
                  levels = c("male", "female"))
)
rownames(sampC) <- colnames(count.mat)

# create rowData
rowData <- DataFrame(
Feature_gr1 = factor(rep(c("blue", "red"), each = 5)),
Feature_gr2 = rep(c(2,3,4,6,1), 2)
)
rownames(rowData) <- rownames(count.mat)

# base tse
tse <- TreeSummarizedExperiment(
  assays = list(counts = count.mat),
  colData = sampC,
  rowData = rowData,
  rowTree = tinyTree
)

# add altExp
alt_colData <- sampC
alt_colData$altTreat <- rep(0:1, 5)

alt_counts <- count.mat[1:5, ]
alt_rowData <- S4Vectors::DataFrame(feature_type = rep("altFeature", 5))
alt_se <- TreeSummarizedExperiment(assays = list(counts = alt_counts), rowData = alt_rowData, colData = alt_colData)

altExp(tse, "altSubset") <- alt_se
metadata(tse)$info <- "test_metadata"
  
  return(tse)
}


# DataFrameSpecs tests

test_that("DataFrameSpecs returns correct structure", {
  tse <- make_test_tse()
  specs <- DataFrameSpecs(colData(tse))
  
  expect_s3_class(specs, "data.frame")
  expect_true(all(c("colName", "colClass", "fctLevels") %in% colnames(specs)))
  
  # condition should be factor with correct levels
  cond_row <- specs[specs$colName == "condition", ]
  expect_equal(cond_row$colClass, "factor")
  expect_true(all(c("trt", "control") %in% strsplit(cond_row$fctLevels, ", ")[[1]]))
  
  # same test but on rowData
  specs <- DataFrameSpecs(rowData(tse))
  
  expect_s3_class(specs, "data.frame")
  expect_true(all(c("colName", "colClass", "fctLevels") %in% colnames(specs)))
  
  # Feature_gr1 should be factor with correct levels
  cond_row <- specs[specs$colName == "Feature_gr1", ]
  expect_equal(cond_row$colClass, "factor")
  expect_true(all(c("blue", "red") %in% strsplit(cond_row$fctLevels, ", ")[[1]]))
})

# test I/O of the _noAltExp() version of function

test_that("TSE round-trips without altExp", {
  tse <- make_test_tse()
  
  tmpDest <- file.path(tempdir(), "testTSE")
  expect_error(write_TSE_to_dir_noAltExp(tse, tmpDest))
  # create directory because the _noAltExp directory does not create it
  dir.create(tmpDest)
  write_TSE_to_dir_noAltExp(tse, tmpDest)
  
  expect_true(file.exists(file.path(tmpDest, "colData.tsv")))
  expect_true(file.exists(file.path(tmpDest, "rowData.tsv")))
  expect_true(file.exists(file.path(tmpDest, "rowTree.tre")))
  expect_false(file.exists(file.path(tmpDest, "colTree.tre")))
  expect_false(dir.exists(file.path(tmpDest, "altExps")))
  
  tse2 <- read_TSE_from_dir_noAltExp(tmpDest)
  
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
  expect_equal(rowLinks(tse2), rowLinks(tse2))
  expect_equal(rowLinks(tse), rowLinks(tse2))
  expect_true(all.equal(rowLinks(tse2), rowLinks(tse2)))
  
  # expect that AltExps are not saved
  expect_true(length(altExpNames(tse)) > 0 & length(altExpNames(tse2)) == 0)
  
  # clean up before next test
  unlink(tmpDest, recursive = TRUE)
})


# test full I/O functions (including AltExp)

test_that("TSE round-trips with altExp", {
  tse <- make_test_tse()
  tmpDest <- file.path(tempdir(), "testTSE")
  write_TSE_to_dir(tse, tmpDest)
  
  expect_true(file.exists(file.path(tmpDest, "colData.tsv")))
  expect_true(file.exists(file.path(tmpDest, "rowData.tsv")))
  expect_true(file.exists(file.path(tmpDest, "rowTree.tre")))
  expect_false(file.exists(file.path(tmpDest, "colTree.tre")))
  # in contrast with the example above, now we expect a TRUE here below
  expect_true(dir.exists(file.path(tmpDest, "altExps")))
  
  tse2 <- read_TSE_from_dir(tmpDest)
  
  # check altExp presence
  expect_true("altSubset" %in% altExpNames(tse2))
  expect_equal(rownames(altExp(tse, "altSubset")), 
               rownames(altExp(tse2, "altSubset")))
  
  # check metadata persisted
  expect_equal(metadata(tse)$info, metadata(tse2)$info)
  # clean up before next test
  unlink(tmpDest, recursive = TRUE)
})

# Expected errors in blank directory rules

test_that("write_TSE_to_dir detects existing directory", {
  tse <- make_test_tse()
  tmpDest <- file.path(tempdir(), "testTSE")
  
  # expect an error on the _noAltExp version of the code
  expect_error(write_TSE_to_dir_noAltExp(tse, tmpDest))
  
  # expect no error on the full code
  expect_no_error(write_TSE_to_dir(tse, tmpDest))
  
  # expect an error if the directory already exists, 
  # even if it's empty
  tmpDest2 <- file.path(tempdir(), "testTSE2")
  expect_error(write_TSE_to_dir(tse, tmpDest))
  
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

