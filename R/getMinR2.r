#' Function that calculates (1 - r-square)
#'
#' This function calculates (1 - r-square) for a set of predicted values (predY) relative to actual empirical values (y).
#' @param predY a vector of numbers that is the fitted data generated from a model.
#' @param y a vector of numbers that were fitted by a model to produce predY.
#' @param na.rm a boolean that specifies whether to remove NAs. DEFAULT = T.
#' @return (1-r2) this is used as the value to be minimized by an optimization program
#' @keywords minimize r2 r-square
#' @export
#' @examples getMinR2 (fittedY, y)

getMinR2 <- function(predY, y, na.rm=T) {

  #if predicted Y is not all NAs, then ...
    if(length(na.omit(predY)) > 1) {
      meanY <- mean(y, na.rm=na.rm)
      ssE <- sum((y - predY)^2, na.rm=na.rm)
      ssT <- sum((y - meanY)^2, na.rm=na.rm)
      r2 <- 1 - (ssE/ssT)
    } else {
      r2 <- 0
    }
      return(1 - r2)
}
