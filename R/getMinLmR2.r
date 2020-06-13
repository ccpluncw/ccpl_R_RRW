#' Function that calculates a linear model and returns (1 - r-square)
#'
#' This function calculates a linear model and returns (1 - r-square).
#' @param x a vector of numbers that is the predictor variable for a linear model.
#' @param y a vector of numbers that the outcome variable for a linear model.
#' @return (1-r2) this is used as the value to be minimized by an optimization program
#' @keywords minimize r2 r-square lm
#' @export
#' @examples min.lm.r2 (x, y)

getMinLmR2 <- function(x, y) {
    #fit a linear model and return (1-r2) if the model fit, otherwise return 1
    out <- 1
    tryCatch ({
      fit.lm <- lm(y~x, data = df.out2)
      out <- (1-summary(fit.lm)$r.squared)
    }, error = function(e) {})

    return(out)
}
