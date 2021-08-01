#' Function that calculates (1 - r-square)
#'
#' This function calculates (1 - r-square) for a set of predicted values (predY) relative to actual empirical values (y).
#' @param predY a vector of numbers that is the fitted data generated from a model.
#' @param y a vector of numbers that were fitted by a model to produce predY.
#' @param standardize A boolean that specifies whether to standardize the scores before calculating the BIC. Standardization is useful if you want to combine different DVs into a single BIC (e.g., RT and accuracy). DEFAULT = FALSE.
#' @return (1-r2) this is used as the value to be minimized by an optimization program
#' @keywords minimize r2 r-square
#' @export
#' @examples getMinR2 (fittedY, y)

getMinR2 <- function(predY, y, standardize = FALSE) {

  #if predicted Y is not all NAs, then ...
    if(length(na.omit(predY)) > 1) {
      r2 <- chutils::ch.R2(y = y, fitY = predY, standardize = standardize)
    } else {
      r2 <- 0
    }
      return(1 - r2)
}
