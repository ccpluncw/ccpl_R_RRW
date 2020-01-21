#' Function that calculates (1 - r-square)
#'
#' This function calculates (1 - r-square) for a set of predicted values (predY) relative to actual empirical values (y).
#' @param predY a vector of numbers that is the fitted data generated from a model.
#' @param y a vector of numbers that were fitted by a model to produce predY.
#' @param na.rm a boolean that specifies whether to remove NAs. DEFAULT = T.
#' @return (1-r2) this is used as the value to be minimized by an optimization program
#' @keywords minimize r2 r-square
#' @export
#' @examples min.r2 (fittedY, y)

min.r2 <- function(predY, y, na.rm=T) {

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

#' Function that calculates a linear model and returns (1 - r-square)
#'
#' This function calculates a linear model and returns (1 - r-square).
#' @param x a vector of numbers that is the predictor variable for a linear model.
#' @param y a vector of numbers that the outcome variable for a linear model.
#' @return (1-r2) this is used as the value to be minimized by an optimization program
#' @keywords minimize r2 r-square lm
#' @export
#' @examples min.lm.r2 (x, y)

min.lm.r2 <- function(x, y) {
    #fit a linear model and return (1-r2) if the model fit, otherwise return 1
    out <- 1
    tryCatch ({
      fit.lm <- lm(y~x, data = df.out2)
      out <- (1-summary(fit.lm)$r.squared)
    }, error = function(e) {})

    return(out)
}
