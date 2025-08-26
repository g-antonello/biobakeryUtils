

#' Extract summary info from `mediation::mediate` model
#'
#' @param mediation.model \code{mediation} object
#'
#' @returns \code{data.frame} with stats as come out of `summary(mediation.model)`
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' #---------- taken from the mediation::mediate first example
#' library(mediation)
#' b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
#' c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)
#' # Estimation via quasi-Bayesian approximation
#' contcont <- mediate(b, c, sims=50, treat="treat", mediator="job_seek")
#' summary(contcont)
#' 
#' # summarized as data.frame
#' tidy_mediation(contcont)
#' }
tidy_mediation <- function(mediation.model){
  med_summary <- summary(mediation.model)
  results_df <- data.frame(
    Effect = c("Average Causal Mediated Effect (ACME)", "Average Direct Effect (ADE)", "Total Effect", "Prop. Mediated"),
    Estimate = c(med_summary$d.avg, med_summary$z.avg, med_summary$tau.coef, med_summary$n.avg),
    CI_Lower_95 = c(med_summary$d.avg.ci[1], med_summary$z.avg.ci[1], med_summary$tau.ci[1], med_summary$n.avg.ci[1]),
    CI_Upper_95 = c(med_summary$d.avg.ci[2], med_summary$z.avg.ci[2], med_summary$tau.ci[2], med_summary$n.avg.ci[2]),
    P_Value = c(med_summary$d.avg.p, med_summary$z.avg.p, med_summary$tau.p, med_summary$n.avg.p)
  )
  
  return(results_df)
}
