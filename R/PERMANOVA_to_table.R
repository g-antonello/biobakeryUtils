#' PERMANOVA to table
#'
#' A utility function that produces a data.frame with PERMANOVA and betadisper
#' statistics from a TreeSummarizedExperiment object
#'
#' @param tse \code{TreeSummarizedExperiment} object that has been already 
#' populated with `mia::addPERMANOVA`
#' @param name \code{character} specifying the `name` paramenter that was passed 
#'to `mia::addPERMANOVA` as name to store the test in the `metadata` slot
#' 
#' @returns Depends on the \code{homogeneity.test} used in
#' \code{mia::addPERMANOVA}:
#' \itemize{
#'   \item no homogeneity test (\code{test.homogeneity = FALSE}): a
#'   \code{data.frame} with the PERMANOVA statistics, one row per term.
#'   \item \code{"permanova"} (default) or \code{"anova"}: the same
#'   \code{data.frame}, with \code{HOMOGENEITY_P},
#'   \code{HOMOGENEITY_TotVariance} and \code{HOMOGENEITY_ExplVariance}
#'   added for each term that was tested for homogeneity (\code{NA} otherwise).
#'   \item \code{"tukeyhsd"}: a \code{list} with \code{permanova} (the
#'   PERMANOVA \code{data.frame}) and \code{homogeneity_pairwise}, a
#'   \code{data.frame} of pairwise group comparisons (\code{comparison},
#'   \code{diff}, \code{lwr}, \code{upr}, \code{p_adj}). Tukey HSD has no
#'   single test per term, so it cannot be merged into the PERMANOVA table.
#'   If several variables were tested, their comparisons are stacked and the
#'   variable is not recorded by \code{mia}.
#' }
#'
#' @importFrom S4Vectors metadata
#' @export 
#'
#' @examples
#' library(mia)
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Apply relative transformation
#' tse <- transformAssay(tse, method = "relabundance")
#' # Perform PERMANOVA
#' tse <- addPERMANOVA(
#'   tse,
#'   assay.type = "relabundance",
#'   method = "bray",  
#'   formula = x ~ SampleType,
#'   permutations = 99
#' )
#' # The results are stored to metadata
#' metadata(tse)[["permanova"]]
#' 
#' # PERMANOVA_to_table makes it faster to retrieve all pieces of it in the
#' # tse object
#' permanova.df <- PERMANOVA_to_table(tse, name = "permanova")
#' permanova.df


PERMANOVA_to_table <- function(tse, name = "permanova") {

    permanova_slot <- metadata(tse)[[name]]

    if (is.null(permanova_slot)) {
      stop("No element named '", name, "' in metadata(tse). ",
           "Was mia::addPERMANOVA run with the same 'name'?", call. = FALSE)
    }

    # when test.homogeneity = FALSE, mia stores the adonis2 table directly,
    # otherwise a list with $permanova and $homogeneity
    if (is.data.frame(permanova_slot)) {
      permanova_tab <- permanova_slot
      homogeneity_tab <- NULL
    } else {
      permanova_tab <- permanova_slot$permanova
      homogeneity_tab <- permanova_slot$homogeneity
    }

    # get PERMANOVA data frame, this is always generated
    permanova_res.df <- data.frame(
      term = rownames(permanova_tab),
      PERMANOVA_df = permanova_tab$Df,
      PERMANOVA_SumOfSqs = permanova_tab$SumOfSqs,
      PERMANOVA_R2 = permanova_tab$R2,
      PERMANOVA_P = permanova_tab[["Pr(>F)"]])

    # case 1 - homogeneity was not tested
    if (is.null(homogeneity_tab)) {
      return(permanova_res.df)
    }

    # case 2 - homogeneity tested with Tukey HSD: pairwise table, one row per
    # pair of groups, which cannot be merged with the per-term PERMANOVA table
    if ("p adj" %in% colnames(homogeneity_tab)) {
      homogeneity_pairwise.df <- data.frame(
        comparison = rownames(homogeneity_tab),
        diff = homogeneity_tab[, "diff"],
        lwr = homogeneity_tab[, "lwr"],
        upr = homogeneity_tab[, "upr"],
        p_adj = homogeneity_tab[, "p adj"],
        row.names = NULL
      )
      return(list(permanova = permanova_res.df,
                  homogeneity_pairwise = homogeneity_pairwise.df))
    }

    # case 3 - homogeneity tested with permutest ("permanova") or anova:
    # one row per variable, p-value in "Pr(>F)" for both
    homogeneity_res.df <- data.frame(
      term = rownames(homogeneity_tab),
      HOMOGENEITY_P = homogeneity_tab[["Pr(>F)"]],
      HOMOGENEITY_TotVariance = homogeneity_tab$`Total variance`,
      HOMOGENEITY_ExplVariance = homogeneity_tab$`Explained variance`
    )

    results_clean <- merge(permanova_res.df, homogeneity_res.df, by = "term", all.x = TRUE)
    # keep row order of permanova results
    results_clean <- results_clean[match(permanova_res.df$term, results_clean$term), ]
    rownames(results_clean) <- NULL

    return(results_clean)
}
