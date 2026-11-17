# shared fixture: small, fast PERMANOVA on GlobalPatterns
make_permanova_tse <- function(...) {
  data("GlobalPatterns", package = "mia", envir = environment())
  tse <- mia::transformAssay(GlobalPatterns, method = "relabundance")
  set.seed(1)
  mia::addPERMANOVA(
    tse,
    assay.type = "relabundance",
    method = "bray",
    formula = x ~ SampleType,
    permutations = 99,
    ...
  )
}

test_that("default homogeneity test ('permanova') reports the permutation p-value", {
  skip_if_not_installed("mia")
  tse <- make_permanova_tse()
  res <- PERMANOVA_to_table(tse)
  homogeneity <- S4Vectors::metadata(tse)$permanova$homogeneity

  expect_s3_class(res, "data.frame")
  expect_identical(res$term[1], "SampleType")

  p <- res$HOMOGENEITY_P[res$term == "SampleType"]
  expect_equal(p, homogeneity["SampleType", "Pr(>F)"])
  expect_true(p > 0 && p <= 1)
  # regression: HOMOGENEITY_P used to be a copy of the explained variance
  expect_false(isTRUE(all.equal(p, res$HOMOGENEITY_ExplVariance[res$term == "SampleType"])))

  # rows without a homogeneity test (Residual, Total) get NA
  expect_true(all(is.na(res$HOMOGENEITY_P[res$term %in% c("Residual", "Total")])))
})

test_that("'anova' homogeneity test reports the F-test p-value", {
  skip_if_not_installed("mia")
  tse <- make_permanova_tse(homogeneity.test = "anova")
  res <- PERMANOVA_to_table(tse)
  homogeneity <- S4Vectors::metadata(tse)$permanova$homogeneity

  expect_equal(res$HOMOGENEITY_P[res$term == "SampleType"],
               homogeneity["SampleType", "Pr(>F)"])
})

test_that("'tukeyhsd' homogeneity test returns PERMANOVA + pairwise tables", {
  skip_if_not_installed("mia")
  tse <- make_permanova_tse(homogeneity.test = "tukeyhsd")
  res <- PERMANOVA_to_table(tse)
  homogeneity <- S4Vectors::metadata(tse)$permanova$homogeneity

  expect_type(res, "list")
  expect_named(res, c("permanova", "homogeneity_pairwise"))
  expect_false(any(grepl("^HOMOGENEITY", colnames(res$permanova))))

  pw <- res$homogeneity_pairwise
  expect_named(pw, c("comparison", "diff", "lwr", "upr", "p_adj"))
  # SampleType has 9 levels -> choose(9, 2) comparisons
  expect_equal(nrow(pw), choose(9, 2))
  expect_equal(pw$p_adj, unname(homogeneity[, "p adj"]))
})

test_that("works when homogeneity was not tested", {
  skip_if_not_installed("mia")
  tse <- make_permanova_tse(test.homogeneity = FALSE)
  res <- PERMANOVA_to_table(tse)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "PERMANOVA_df", "PERMANOVA_SumOfSqs", "PERMANOVA_R2", "PERMANOVA_P"))
})

test_that("informative error when the metadata name does not exist", {
  skip_if_not_installed("mia")
  tse <- make_permanova_tse()
  expect_error(PERMANOVA_to_table(tse, name = "not_there"), "addPERMANOVA")
})
