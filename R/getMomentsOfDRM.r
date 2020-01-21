#' This function fits a  3 parameter log-Logistic to a Cumulative Distribution Function, and returns the moments of the function.
#'
#' This function fits a  3 parameter log-Logistic to a Cumulative Distribution Function, and returns the moments of the function.
#' @param CDFps A vector containing the proportions of the CDF.
#' @param values A vector containing the values (x-axis) of the CDF.
#' @return A dataframe that contains the moments of the simmulated values from the CDF ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the probability of crossing the boundary)
#' @keywords DRM CDF moments values
#' @export
#' @examples getMomentsOfDRM (CDFvector, valueVector)

getMomentsOfDRM <- function(CDFps, values) {

  tmp.drm <- NULL
  df.out <- NULL

  if( mean(CDFps[trunc((length(CDFps) - length(CDFps)/10)):length(CDFps)]) > 0.005) {
    tryCatch ({
      tmp.drm <- drc::drm(CDFps ~ values, fct = LL2.3(), robust='mean')
      }, error = function(e) {
#          print(paste("drm function did not fit - trying alternate method ..."))
    })
  }
  if (is.null(tmp.drm)) {
    pdf <- diff(CDFps)/diff(values)
    pdf <- ifelse(pdf < 0, 0, pdf)
    df.out <- NULL
    tryCatch ({
      df.out <- getMomentsOfPDF(pdf, values[2:length(values)])
      }, error = function(e) {
#          print(paste("getMomentsOfPDF - trying alternate method ..."))
    })
    if(is.null(df.out)) {
      df.out$mean <- NA
      df.out$sd <- NA
      df.out$Q25 <- NA
      df.out$Q50 <- NA
      df.out$Q75 <- NA
    }
    df.out$pCross <- mean(CDFps[(length(CDFps) - length(CDFps)/10):length(CDFps)])
    if(df.out$pCross < 0.01) {
      df.out$mean <- NA
      df.out$sd <- NA
      df.out$Q25 <- NA
      df.out$Q50 <- NA
      df.out$Q75 <- NA
    }
  } else {
      pdf <- diff(predict(tmp.drm))/diff(values)
      pdf <- ifelse(pdf < 0, 0, pdf)
      df.out <- NULL
      tryCatch ({
        df.out <- getMomentsOfPDF(pdf, values[2:length(values)])
        }, error = function(e) {
#            print(paste("getMomentsOfPDF - trying alternate method ..."))
      })
      if(is.null(df.out)) {
        df.out$mean <- NA
        df.out$sd <- NA
        df.out$Q25 <- NA
        df.out$Q50 <- NA
        df.out$Q75 <- NA
      }
      df.out$pCross <- max(fitted(tmp.drm))
  }

  return(df.out)

}
