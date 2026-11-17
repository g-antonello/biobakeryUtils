#' Single-outcome Meta-analysis with metafor
#' 
#' @description
#' This function is a simple wrapper for what I need of a typical meta-analysis
#' workflow for micobiome data. This function is called by \code{metaAnalyze},
#' which generalizes its use.
#' 
#'
#' @param yi \code{numeric} vector of effect sizes `beta`
#' @param sei \code{numeric} vector of standard errors of `beta`
#' @param method \code{chracter}. One of the options from the `method` parameter
#'   in metafor::rma. (Default: 'REML')
#' @param test \code{chracter}. One of the options from the `test` parameter
#'   in metafor::rma. (Default: 'knha', which is suggested for meta-analyses 
#'   with up to 20 datasets.
#'
#' @returns \code{numeric} vector containing the meta-analysis basic stats:
#' StudiesTested: number of studies with feature, Beta.pool: random effects beta
#' SE.pool: standard error of the Beta.pool point esimate; CI.L: 95% lower bound 
#' confidence interval of the Beta and SE; CI.R: upper 95% bound confidence 
#' interval; P.Value: meta-analysis p-value estimated with the method in `test`;
#' Het.tau2: tau-squared heterogeneity estimate, Het.I2: I-squared heterogeneity
#' estimate Het.Pval: P-value associated with heterogeneity.
#' @export
#' 
#' @importFrom metafor rma
#' 
#' @examples
#' library(metafor)
#' dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#' metaAnalyzeSingle(yi = dat$yi, sei = dat$vi, method = "REML")
metaAnalyzeSingle <- function(yi, sei, method = "REML", test = "knha") {
  
  # make fit. the REML method is the gold standard for random effects meta-analyses
  # the "knha" (Knapp-Hartung adjustment) instead of "z" because z requires at
  # least 20 cohorts
  if(is.null(test)){
    test <- ifelse(length(yi) >= 20, "z", "knha")
  }
  
  fit <- tryCatch(
    suppressWarnings(metafor::rma(yi = yi, sei = sei, method = method, test = test)),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  return(
    c(
      StudiesTested = fit$k,
      Beta.pool = fit$beta[[1]],
      SE.pool   = fit$se,
      CI.L      = fit$ci.lb,
      CI.R      = fit$ci.ub,
      P.Value   = fit$pval,
      Het.tau2      = fit$tau2,
      Het.I2        = fit$I2,
      Het.Qstat         = fit$QE,
      Het.Pval    = fit$QEp
    )
  )
}

#---------------------------------------------------------------------------

#' Meta-Analyze multi-group 
#'
#' @param res.df \code{data.frame} with long-format meta-analysis results 
#' @param grouping.col \code{character} variable that was used to group data
#'   which were later meta-analyzed
#' @param beta.col \code{character}. Column name of the `beta` estimates
#' @param se.col \code{character}. Column name of the `SE` estimates
#' @param method \code{character}. Meta-analysis method as found in the omonimous
#'   parameter in metafor::rma
#' @param test  \code{character}. Meta-analysis test as found in the omonimous
#'   parameter in metafor::rma
#' @param minStudies.n \code{integer}. Minimum number of studies to retain a 
#'   feature for meta-analysis 
#' @param p.adj_method \code{character}. One of the possible options in 
#'   `p.adjust`
#'
#' @returns A \code{list} of 3 elements: (1) "Beta_matrix"; (2) "SE_matrix"; 
#'  (3) "MetaAnalysis". The first 2 are the matrices used to calculate the meta-
#'  analyses, while the 3rd one is a `data.frame` of all results together
#' @export
#' 
#' @importFrom reshape2 acast
#' @importFrom stats p.adjust as.formula
#' @examples
#' library(curatedMetagenomicData)
#' library(mia)
#' library(magrittr)
#' library(dplyr)
#' 
#' sex_meta_analysis <- sampleMetadata %>%
#'   # keep cases where we know sex and that are not longitudinal samples
#'   filter(
#'     !is.na(sampleMetadata$gender),
#'     is.na(sampleMetadata$visit_number),
#'     body_site == "stool", 
#'     disease == "healthy"
#'   )
#' 
#' # restrict to datasets with at least 30 males and 30 females
#' ok_dataset <- sex_meta_analysis %>%
#'   group_by(study_name) %>%
#'   reframe(
#'     nMale = sum(gender == "male"),
#'     nFemale = sum(gender == "female")
#'   ) %>%
#'   filter(nMale > 30, nFemale > 30)
#' 
#' sex_meta_analysis <- filter(sex_meta_analysis, study_name %in% ok_dataset$study_name)
#' sex_meta_analysis.se <- returnSamples(sex_meta_analysis, 
#'                                       dataType = "relative_abundance", counts = TRUE)
#' 
#' # make a copy of the assay, it actually has counts, as requested
#' # testable with `range(assay(sex_meta_analysis.se))`
#' assay(sex_meta_analysis.se, "counts") <- assay(sex_meta_analysis.se, 
#'                                                "relative_abundance")
#' 
#' 
#' # filter by prevalence, at least 5 counts in each cell, at least 1% of 
#' # the samples
#' sex_meta_analysis.se <- mia::subsetByPrevalent(sex_meta_analysis.se, 
#'   assay.type = "counts", detection = 5, prevalence = 0.01)
#' # fix error of same feature in multiple lines by aggregating
#' 
#' sex_meta_analysis.se <- mia::agglomerateByRank(sex_meta_analysis.se, 
#'   "species")
#' 
#' 
#' singledataset_analysis_sex.df <- sex_meta_analysis.se %>%
#'   splitOn(group = "study_name") %>%
#'   lapply(function(x.tse) {
#'     x.tse <- x.tse[rowSums(assay(x.tse, "relative_abundance")) > 0,]
#'     return(limmaTSE(tse = x.tse, assay.type = "counts", 
#'                       voom.transform = TRUE, 
#'                       formula = ~ gender, coef = "gendermale", 
#'                       p.adj_method = "BH"))
#'   }
#'   ) %>% bind_rows(.id = "study_name")
#' 
#' meta_analysis_results <- metaAnalyze(res.df = singledataset_analysis_sex.df,
#'                                      grouping.col = "study_name", 
#'                                      beta.col = "Beta", 
#'                                      se.col = "SE", method = "REML", 
#'                                      test = "knha", 
#'                                      minStudies.n = 6, 
#'                                      p.adj_method = "BH")
#' 
#' str(meta_analysis_results)

metaAnalyze <- function(res.df,
                        grouping.col = "study_name",
                        beta.col = "Beta",
                        se.col = "SE",
                        method = "REML",
                        test = "knha",
                        minStudies.n = 6,
                        p.adj_method = "BH") {
  
  avail_methods <- c("DL", "HE", "CO", "HS", "HSk", "SJ", "ML", "REML", "EB", "PM", "GENQ", "PMM", "GENQM", "FE")
  avail_tests <- c("z", "knha", "t")
  FeatureID_indx <- which(colnames(res.df) == "N") - 1
  FeatureID <- colnames(res.df)[FeatureID_indx]
  
  # rest of rowData information, only unique rows
  restOfrowData <- res.df[,2:FeatureID_indx]
  restOfrowData <- restOfrowData[!duplicated(restOfrowData),]
  
  if(!(method %in% avail_methods)){
    stop("`method` requested not available. Choose between any found under 
         'method' in the metafor::rma documentation. E.g.: ", 
         paste(avail_methods, collapse = "; "))
  }
  
  if(!(test %in% avail_tests)){
    stop("`test` requested not available. Choose between: ", paste(avail_tests, 
                                                                   collapse = "; "))
  }
  
  beta.mtx <- reshape2::acast(res.df, formula = as.formula(sprintf("%s ~ %s", grouping.col, FeatureID)), value.var = beta.col)
  se.mtx <- reshape2::acast(res.df, formula = as.formula(sprintf("%s ~ %s", grouping.col, FeatureID)), value.var = se.col)
  
  # Drop features that were tested in less than the minStudies.n wanted
  beta.mtx_ok <- beta.mtx[, colSums(is.na(beta.mtx)) < minStudies.n]
  se.mtx_ok <- se.mtx[, colSums(is.na(se.mtx)) < minStudies.n]
  
  results.list <- lapply(seq_len(ncol(beta.mtx_ok)), function(i)
    metaAnalyzeSingle(
      yi = beta.mtx_ok[, i],
      sei = se.mtx_ok[, i],
      method = method,
      test = test
    ))
  names(results.list) <- colnames(beta.mtx_ok)
  
  # format results as data.frame
  results_formatted <- as.data.frame(do.call(rbind, results.list))
  results_formatted <- cbind(rownames(results_formatted), results_formatted)
  colnames(results_formatted)[1] <- FeatureID
  rownames(results_formatted) <- NULL
  
  # adjust p-values and reorder columns
  results_formatted$P.Value.adj = p.adjust(results_formatted$P.Value, method = p.adj_method)
  results_formatted <- results_formatted[, c(1:7, 12, 8:11)]
  
  # merge rest of rowData from the meta-analysis and preserve column order
  results_formatted_with_rowData <- merge(restOfrowData, results_formatted, by = FeatureID)
  results_formatted_with_rowData <- results_formatted_with_rowData[, c(colnames(restOfrowData), setdiff(colnames(results_formatted), colnames(restOfrowData)))]
  # reorder by meta-analysis adjusted p-value
  results_formatted_with_rowData <- results_formatted_with_rowData[order(results_formatted_with_rowData$P.Value.adj), ]
  
  # make a beta matrix with study name as first column
  beta.mtx_ok <- cbind.data.frame(rownames(beta.mtx_ok), beta.mtx_ok)
  colnames(beta.mtx_ok)[1] <- "Dataset"
  # make a se matrix with study name as first column
  se.mtx_ok <- cbind.data.frame(rownames(se.mtx_ok), se.mtx_ok)
  colnames(se.mtx_ok)[1] <- "Dataset"
  
  return(
    list(
      "Beta_matrix" = beta.mtx_ok,
      "SE_matrix" = se.mtx_ok,
      "MetaAnalysis" = results_formatted_with_rowData
    )
  )
}
