#' Function that calculates AIC for the entire RRW dataset
#'
#' This function calculates AIC for a set of predicted values (predY) relative to actual empirical values (y).
#' @param simPhit a vector of numbers that is the fitted data generated from the model for the probability correct.
#' @param dataPhit a vector of numbers that were fitted by a model to produce simPhit.
#' @param simRT a vector of numbers that is the fitted data generated from the model for the RT.
#' @param dataPhit a vector of numbers that were fitted by a model to produce simRT.
#' @param numParameters The number of free parameters.
#' @param standardize AIC is a function of the number of observations. If you want to correct for the number of observations, then set "standardize" to TRUE. DEFAULT = FALSE.
#' @return AIC this is used as the value to be minimized by an optimization program
#' @keywords minimize AIC
#' @export
#' @examples getAICRRW (fittedPhit, pHit, fittedRT, rt, 5, standardize = TRUE)

getAICRRW <- function (simPhit, dataPhit, simRT, dataRT, numParameters, standardize = FALSE) {
  sim <- c(scale(simPhit), scale(simRT))
  dat <- c(scale(dataPhit), scale(dataRT))

  out.AIC.final <- chutils::ch.AIC(dat, sim, numParameters, standardize = standardize)

  return(out.AIC.final)
}
