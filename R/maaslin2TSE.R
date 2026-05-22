#' Maaslin2 Wrapper for TreeSummarizedExperiment (TSE) Objects
#'
#' @description 
#' This function provides a seamless wrapper to run MaAsLin2 
#' directly on a `TreeSummarizedExperiment` object. It also includes the benefit
#' of sinking console output while running it, and deletes all temporary files 
#' if requested (default 'yes'). Finally, the function return a cleaned, 
#' formatted data frame containing the statistical results merged with 
#' feature metadata and some per-feature statistics (rowData).
#' 
#' @param tse A `TreeSummarizedExperiment` or `SummarizedExperiment` object 
#'   containing biobakery-derived tabular data (MetaPhlAn, HUMAnN).
#' @param assay.type Character string specifying the name of the assay to use 
#'   (default: "relative_abundance").
#' @param outdir Character string specifying the output directory for temporary 
#'   MaAsLin2 files. Files are unlinked after execution (default: "output/maaslin2").
#' @param outdir_delete \code{logical}. Should the output directory be deleted?
#'   default = TRUE
#' @param min_abundance Numeric. The minimum abundance threshold for filtering 
#'   features (default: 0.0).
#' @param min_prevalence Numeric. The minimum prevalence threshold for filtering 
#'   features. (default: 0.1).
#' @param min_variance Numeric. The minimum variance threshold for filtering 
#'   features (default: 0.0).
#' @param normalization Character. The normalization method to apply. 
#'   Options: "TSS", "CLR", "CSS", "NONE", "TMM" (default: "TSS").
#' @param transform Character. The transformation method to apply. 
#'   Options: "LOG", "LOGIT", "AST", "NONE" (default: "LOG").
#' @param analysis_method Character. The analytical method to use. 
#'   Options: "LM", "CPLM", "ZINB", "NEGBIN" (default: "LM").
#' @param max_significance Numeric. The q-value threshold for significance 
#'   (default: 0.25).
#' @param random_effects Character vector indicating the metadata columns to be 
#'   used as random effects (default: NULL).
#' @param fixed_effects Character vector indicating the metadata columns to be 
#'   used as fixed effects. Must be present in `colData(tse)`.
#' @param correction Character. Multiple hypothesis testing correction method. 
#'   Options include "BH", "holm", etc. (default: "BH").
#' @param standardize Logical. Whether to standardize continuous metadata 
#'   (default: TRUE).
#' @param cores Integer. Number of cores for parallel processing (default: 1).
#' @param plot_heatmap Logical. Whether to generate a heatmap plot (default: TRUE).
#' @param heatmap_first_n Integer. Number of top features to include in the heatmap 
#'   (default: 50).
#' @param plot_scatter Logical. Whether to generate scatter plots for significant 
#'   associations (default: TRUE).
#' @param max_pngs Integer. Maximum number of PNGs to save (default: 10).
#' @param save_scatter Logical. Whether to save scatter plots to disk (default: FALSE).
#' @param save_models Logical. Whether to save the fitted models (default: FALSE).
#' @param reference Character string or vector specifying the reference level(s) 
#'   for categorical fixed effects (default: NULL).
#'
#' @returns A data frame containing the merged results of MaAsLin2 statistics 
#'   (Beta, SE, t_stat, Confidence Intervals, P-values, and adjusted P-values) 
#'   alongside the feature metadata (`rowData`) and compositional statistics. 
#'   Results are ordered by adjusted p-value.
#' 
#' @importFrom utils read.delim capture.output
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom Maaslin2 Maaslin2
#' @export
#'
#' @examples
#' library(mia)
#' data(enterotype)
#' enterotype <- transformAssay(enterotype, 
#'   assay.type = "counts", 
#'   method = "clr", 
#'   pseudocount = TRUE)
#' 
#' no_NA_variables <- complete.cases(colData(enterotype)[,c("Enterotype","Age")])
#' enterotype_noNAs <- enterotype[,no_NA_variables]
#' DA_results <- maaslin2TSE(enterotype_noNAs, 
#'   assay.type = "counts",
#'   transform = "LOG",
#'   fixed_effects = c("Enterotype","Age"))
#' 
#' str(DA_results)
maaslin2TSE <- function(tse, 
                        assay.type = "relative_abundance", 
                        outdir = NULL, 
                        outdir_delete = TRUE,
                        min_abundance = 0.0,
                        min_prevalence = 0.1,
                        min_variance = 0.0,
                        normalization = "TSS",
                        transform = "LOG", 
                        analysis_method = "LM",
                        max_significance = 0.25,
                        random_effects = NULL,
                        fixed_effects = NULL,
                        correction = "BH",
                        standardize = TRUE,
                        cores = 1,
                        plot_heatmap = TRUE, 
                        heatmap_first_n = 50,
                        plot_scatter = TRUE,
                        max_pngs = 10,
                        save_scatter = FALSE,
                        save_models = FALSE,
                        reference = NULL){
  
  # create outdir 
  if(is.null(outdir)){
    outdir <- tempdir()
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  
  # 2. Extract and format feature data (Samples as rows, Features as columns)
  input_data.df <- as.data.frame(t(SummarizedExperiment::assay(tse, assay.type)))
  
  # 3. Extract metadata, again with samples as rows
  input_metadata.df <- as.data.frame(SummarizedExperiment::colData(tse))
  
  # 4. Handle fixed and random effects defaults
  if (is.null(fixed_effects)) {
    stop("You must specify at least one column from colData(tse) as 'fixed_effects'.")
  }
  
  # 5. Execute Maaslin2 safely by sinking output to a file
    
      Maaslin2::Maaslin2(
      input_data       = input_data.df,
      input_metadata   = input_metadata.df,
      output           = outdir,
      fixed_effects    = fixed_effects,
      random_effects   = random_effects,
      plot_heatmap     = plot_heatmap,
      plot_scatter     = plot_scatter,
      reference        = reference,
      min_abundance    = min_abundance,
      min_prevalence   = min_prevalence,
      min_variance     = min_variance,
      normalization    = normalization,
      transform        = transform,
      heatmap_first_n  = heatmap_first_n,
      analysis_method  = analysis_method,
      max_significance = max_significance,
      correction       = correction,
      standardize      = standardize,
      cores            = cores,
      max_pngs         = max_pngs,
      save_scatter     = save_scatter,
      save_models      = save_models
    )
    
  # 6. Load results and unlink everything else
  results_raw <- read.delim(file.path(outdir, "all_results.tsv"))
  
  # 7. Format results better. getFeatureStats is exported in biobakeryUtils
  compositionalStats <- getFeatureStats(tse, assay.type)
  compositionalStats$feature <- rownames(compositionalStats)
  
  # 7.1 - get rowData and relative compositional stats
  rowData.df <- as.data.frame(SummarizedExperiment::rowData(tse))
  rowData.df$feature <- rownames(rowData.df)
  rowData_with_stats.df <- merge(rowData.df, compositionalStats, by = "feature")
  
  # 7.2 - merge rowData with results and fix columns with same name
  results_with_rowadata <- merge(rowData_with_stats.df, results_raw, by = "feature")
  results_with_rowadata$N.y <- NULL
  colnames(results_with_rowadata)[which(colnames(results_with_rowadata) == "N.x")] <- "N" 
  
  # 7.3 - filter results to include only metadata of interest
  results_with_rowadata_filt <- results_with_rowadata[results_with_rowadata$metadata %in% fixed_effects,]
  
  # 8 - Build final dataframe using the FILTERED object, not the raw one
  results_clean <- data.frame(
    results_with_rowadata_filt[, colnames(rowData_with_stats.df)],
    ExposureValue = paste0(results_with_rowadata_filt$metadata, results_with_rowadata_filt$value),
    Beta        = results_with_rowadata_filt$coef,
    SE          = results_with_rowadata_filt$stderr,
    t_stat      = results_with_rowadata_filt$coef / results_with_rowadata_filt$stderr,
    CI.L        = results_with_rowadata_filt$coef - (results_with_rowadata_filt$stderr * 1.96),
    CI.R        = results_with_rowadata_filt$coef + (results_with_rowadata_filt$stderr * 1.96),
    P.Value     = results_with_rowadata_filt$pval,
    P.Value.adj = results_with_rowadata_filt$qval
  )
  
  results_clean$feature <- NULL
  results_clean <- results_clean[order(results_clean$P.Value.adj),]
  
  # second to last step, delete directory if requested
  if(outdir_delete){
    unlink(outdir, recursive = TRUE)
  }
  
  # 9. Return the clean results data frame
  return(results_clean)
}
