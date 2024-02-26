#' Function that calculates (1 - r-square) for the entire RRW dataset
#'
#' This function calculates (1 - r-square) for a set of predicted values (predY) relative to actual empirical values (y).
#' @param simPhit a vector of numbers that is the fitted data generated from the model for the probability correct.
#' @param dataPhit a vector of numbers that were fitted by a model to produce simPhit.
#' @param simRT a vector of numbers that is the fitted data generated from the model for the RT.
#' @param dataPhit a vector of numbers that were fitted by a model to produce simRT.
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the r square of the pHit and RT is calculated by equalizing the influence of the number of observations. If it is set to FALSE, then the output is calculated as usual. DEFAULT = FALSE.
#' @param standardize a boolean that specifies whether to standardize the DVs before running. Unstandardized data may give biased results when equalizeRTandPhit = T. DEFAULT = T.
#' @return (1-r2) this is used as the value to be minimized by an optimization program
#' @keywords minimize r2 r-square
#' @export
#' @examples getMinR2RRW (fittedPhit, pHit, fittedRT, rt)

getMinR2RRW <- function(simPhit, dataPhit, simRT, dataRT, equalizeRTandPhit = FALSE, standardize=T) {

  #### put both dv's on the same scale
  if (standardize) {
    df.rt.z <- standardizeDataAndFit(dataRT, simRT)
    dataRT <- df.rt.z$data
    simRT <- df.rt.z$fit
    df.pHit.z <- standardizeDataAndFit(dataPhit, simPhit)
    dataPhit <- df.pHit.z$data
    simPhit <- df.pHit.z$fit
  }
  #####

  #get rss and tss for pHit and rt
  pHit.tss <- chutils::ch.TSS(dataPhit, standardize = FALSE)
  pHit.rss <- chutils::ch.RSS(dataPhit, simPhit, standardize = FALSE)
  if(is.na(pHit.rss)) {
    #if pHit.rss == NA, then set rss equal to tss so r2 = 0
    pHit.rss <- pHit.tss
  }

  rt.tss <- chutils::ch.TSS(dataRT, standardize = FALSE)
  rt.rss <- chutils::ch.RSS(dataRT, simRT, standardize = FALSE)
  if(is.na(rt.rss)) {
    #if rt.rss == NA, then set rss equal to tss so r2 = 0
    rt.rss <- rt.tss
  }

  if (equalizeRTandPhit) {
    #remove the influence of the number of observations from the pHit rss and tss by dividing by n.
    #get length without NAs
    pHit.n <- length(dataPhit[!is.na(dataPhit)])
    pHit.ave.rss <- pHit.rss/pHit.n
    pHit.ave.tss <- pHit.tss/pHit.n

    #and the rt rss
    #get length without NAs
    rt.n <- length(dataRT[!is.na(dataRT)])
    rt.ave.rss <- rt.rss/rt.n
    rt.ave.tss <- rt.tss/rt.n

    total.n <- pHit.n + rt.n

    #Then average the equalized rss from RT and pHit t
    #and, reintroduce the influence of n by multiplying the equalized rss by the total n.
    rss.equal <- mean( c(pHit.ave.rss, rt.ave.rss), na.rm = T ) * total.n
    tss.equal <- mean( c(pHit.ave.tss, rt.ave.tss), na.rm = T ) * total.n

    r2 <- 1 - (rss.equal/tss.equal)

  } else {
    df.rt.z <- standardizeDataAndFit(dataRT, simRT)
    df.pHit.z <- standardizeDataAndFit(dataPhit, simPhit)
    df.z <- rbind(df.rt.z, df.pHit.z)
    r2 <- ch.R2(df.z$data, df.z$fit)
  }
  #return 1-r2 for the minimization.
  out.r2 <- 1-r2

  return(out.r2)
}
