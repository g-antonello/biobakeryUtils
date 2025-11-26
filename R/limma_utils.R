#' limma wrapper for biobakery data in TreeSummarizedExperiment format
#'
#' NB: This particular wrapper does not perforom limma::voom, but simply the 
#' `lmFit` function with a given formula. It does not perform transformations
#' either. IMPORTANT: taxa stats are calculated on the first assay, which is 
#' expected to be a `abundance/relabundance` as clean as possible from either
#' MetaPhlAn or HUMAnN.
#' 
#' @param tse \code{(Tree)SummarizedExperiment}. Input data with available rowData
#' @param assay.type \code{character}. The assay to test with limma
#' @param formula \code{formula}. Model to pass to `model.matrix` before lmFit
#' @param coef \code{character}. the coefficient as seen in one of `limma$coef` 
#' columns. This is not so generalizable yet
#'
#' @returns a \code{data.frame} with a lot of feature information along with 
#' limma's summary statistics
#' 
#' @export
#'
#' @examples
#' data(enterotype)
#' enterotype <- transformAssay(enterotype, assay.type = "counts", method = "clr", pseudocount = TRUE)
#' enterotype_noNAs <- enterotype[,!is.na(enterotype@colData$Enterotype) & !is.na(enterotype@colData$Age)]
#' limmaResults <- limmaFit_TSE(enterotype_noNAs, assay.type = "clr", formula = ~ Enterotype + Age, coef = "Enterotype3")
#' 
#' str(limmaResults)

limmaFit_TSE <- function(tse, assay.type, formula, coef){
  if(is.null(assay.type)){
    assay.type = assayNames(tse)[1]
  }
  
  if(is.null(coef)){
    stop("coef not provided but necessary for table extraction")
  }
  
  # Calculate preliminary stats and colmns to add to output later
  rowData.df <- as.data.frame(rowData(tse))
  rowData.df$feature <- rownames(rowData.df)
  
  compositionalStats <- getFeatureStats(tse, assayNames(tse)[1])
  compositionalStats$feature <- rownames(compositionalStats)
  
  rowData_with_Stats.df <- merge(rowData.df, compositionalStats, by = "feature")
  
  # calculate limma model
  
  modMtx <- model.matrix(formula, data = colData(tse))
  
  lmFitBasic <- lmFit(assay(tse, assay.type), design = modMtx)
  lmFitBasic <- eBayes(lmFitBasic)
  
  lmFit_table <- topTable(lmFitBasic, coef = coef, number = Inf, confint = TRUE, adjust.method = "BH")
  lmFit_table$feature <- rownames(lmFit_table)
  
  final_df <- merge(rowData_with_Stats.df, lmFit_table, by = "feature")
  final_df$SE <- (final_df$CI.R - final_df$CI.L)/(1.96*2)
  # Reorganize columns
  final_df <- final_df[, c(colnames(rowData_with_Stats.df), "logFC","SE", "t", "B", "CI.L", "CI.R", "P.Value", "adj.P.Val")]
  colnames(final_df) <- c(colnames(rowData_with_Stats.df), "Beta","SE", "t_stat", "B_stat", "CI.L", "CI.R", "P.Value", "adj.P.Value")
  
  # remove merging column
  final_df$feature <- NULL
  
  return(final_df)
}

