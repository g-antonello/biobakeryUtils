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
  fields <- list(
    list(label = "Average Causal Mediated Effect (ACME)", est = "d.avg", ci = "d.avg.ci", p = "d.avg.p"),
    list(label = "Average Direct Effect (ADE)",           est = "z.avg", ci = "z.avg.ci", p = "z.avg.p"),
    list(label = "Total Effect",                          est = "tau.coef", ci = "tau.ci", p = "tau.p"),
    list(label = "Prop. Mediated",                        est = "n.avg", ci = "n.avg.ci", p = "n.avg.p")
  )
  
  return(
    data.frame(
    Effect      = vapply(fields, `[[`, character(1), "label"),
    Estimate    = vapply(fields, \(f) mediation.model[[f$est]],    numeric(1)),
    CI_Lower_95 = vapply(fields, \(f) mediation.model[[f$ci]][1], numeric(1)),
    CI_Upper_95 = vapply(fields, \(f) mediation.model[[f$ci]][2], numeric(1)),
    P.value     = vapply(fields, \(f) mediation.model[[f$p]],      numeric(1))
  )
  )
}
