#' This function gets the predicted value of the three parameter log logisistic function returned by drc::drm().
#'
#' This function gets the predicted value of the three parameter log logisistic function returned by drc::drm().
#'
#'
#' @param val A number that is the "x" value for which you want the "y" value returned.
#' @param drmModel The DRM model returned by drc::drm().
#'
#' @return The "Y" value that corresponds to the "x" value of the "val" parameter
#' @keywords DRM log logistic
#' @export
#' @examples getDRMpredValue (drmVal, drmModel)

getDRMpredValue <- function(val, drmModel) {
  b <- as.numeric(drmModel$coefficients["b:(Intercept)"])
  d <- as.numeric(drmModel$coefficients["d:(Intercept)"])
  e <- as.numeric(drmModel$coefficients["e:(Intercept)"])

  out <- d/(1+exp(b*(log(val) - log(e))))
  return(1-out)

}
