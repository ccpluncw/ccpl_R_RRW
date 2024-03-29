#' Function that calculates AIC or BIC for the entire RRW dataset
#'
#' This function calculates AIC or BIC for a set of predicted values (predY) relative to actual empirical values (y).
#' @param simPhit a vector of numbers that is the fitted data generated from the model for the probability correct.
#' @param dataPhit a vector of numbers that were fitted by a model to produce simPhit.
#' @param simRT a vector of numbers that is the fitted data generated from the model for the RT.
#' @param dataPhit a vector of numbers that were fitted by a model to produce simRT.
#' @param numParameters The number of free parameters.
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations in the BIC.  RT has more observations that pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the rss is standardized by the number of observations. If it is set to FALSE, then the BIC is calculated as usual. DEFAULT = FALSE.
#' @param ICtype A string specifying whether to return the AIC or BIC.  Valid inputs are: "AIC" and "BIC".  DEFAULT = "BIC"
#' @param standardize a boolean that specifies whether to standardize the DVs before running. Unstandardized data may give biased results when equalizeRTandPhit = T. DEFAULT = T.
#''
#' @return BIC this is used as the value to be minimized by an optimization program
#' @keywords minimize BIC
#' @export
#' @examples getICRRW (fittedPhit, pHit, fittedRT, rt, 5, standardize = TRUE, ICtype = "BIC")

getICRRW <- function (simPhit, dataPhit, simRT, dataRT, numParameters, equalizeRTandPhit = FALSE, ICtype = "BIC", standardize=T) {

  #make sure minimizeStat is valid
  validOpts <- c("BIC", "AIC")
  if(!(ICtype %in% validOpts) ) {
    stop (paste("you set ICtype to:", ICtype, ", but it must be one of the following:", validOpts, sep=" "))
  }

    #To equalize the influence of pHit fit and RT fit:
  if (equalizeRTandPhit) {
    if (standardize) {
      df.rt.z <- standardizeDataAndFit(dataRT, simRT)
      dataRT <- df.rt.z$data
      simRT <- df.rt.z$fit
      df.pHit.z <- standardizeDataAndFit(dataPhit, simPhit)
      dataPhit <- df.pHit.z$data
      simPhit <- df.pHit.z$fit
    }

    #First remove the influence of the number of observations from the pHit rss by dividing by n.
    pHit.rss <- chutils::ch.RSS(dataPhit, simPhit, standardize = FALSE)
    #if pHit function is not fit, then IC should be infinite
    if(is.na(pHit.rss)) {
      pHit.rss <- Inf
    }
      #get length without NAs
    pHit.n <- length(dataPhit[!is.na(dataPhit)])
    pHit.ave.rss <- pHit.rss/pHit.n

      #and the rt rss
    rt.rss <- chutils::ch.RSS(dataRT, simRT, standardize = FALSE)
    #if RT function is not fit, then IC should be infinite
    if(is.na(rt.rss)) {
      rt.rss <- Inf
    }
      #get length without NAs
    rt.n <- length(dataRT[!is.na(dataRT)])
    rt.ave.rss <- rt.rss/rt.n

    total.n <- pHit.n + rt.n

      #Then average the equalized rss from RT and pHit t
      #and, reintroduce the influence of n by multiplying the equalized rss by the total n.
    rss.equal <- mean( c(pHit.ave.rss, rt.ave.rss), na.rm = T ) * total.n

    out.IC.final <- chutils::ch.ICfromRSS(rss = rss.equal, n = total.n, numParameters = numParameters, ICtype = ICtype)

  } else {
    df.rt.z <- standardizeDataAndFit(dataRT, simRT)
    df.pHit.z <- standardizeDataAndFit(dataPhit, simPhit)
    df.z <- rbind(df.rt.z, df.pHit.z)

    out.IC.final <- chutils::ch.IC(df.z$data, df.z$fit, numParameters = numParameters, standardize = FALSE, ICtype = ICtype)
  }

  return(out.IC.final)
}
