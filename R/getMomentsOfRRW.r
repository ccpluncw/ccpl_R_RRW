#' This function simulates the RRW for a single set of parameters for a single overlap.
#'
#' This function simulates the RRW for a single set of parameters for a single overlap.
#'
#'
#' @param overlap A number between 0-1 representing a single distributional overlap.
#' @param boundary A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param correctBoundary Either 1 or -1.  If correctBoundary = 1, then the positive boundary is the boundary for a correct response.  If correctBoundary = -1, then the negative boundary is the boundary for a correct response.
#' @param stepSize A number that specifies the granularity that the simulation will generate the probability crossing a boundary.  stepSize = 1 is the Default and most granular.
#' @param startValue A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  startValue = 0 is the default and represents an unbiased start point.
#' @param noiseSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. noiseSD = 0 is the default and represents no noise being added to each step.
#' @param decayAsymptote A positive number representing the asymptote value in the Information Accrual Bias (IAB). decayAsymptote = 0.2 is the default.
#' @param decayBeta A signed number representing the beta value in the Information Accrual Bias (IAB). decayBeta = 0 is the default and represents no IAB.
#' @param loops A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#'
#' @return A dataframe that contains the moments of the simmulated values from the model ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the probability of crossing each boundary); correct (TRUE or FALSE)), along with the parameters of the model.
#' @keywords RRW random walk simulation moments overlap
#' @export
#' @examples getMomentsOfRRW (overlap = 0.1, boundary=14, startValue=0.1, loops = 400)


getMomentsOfRRW <- function (overlap, boundary, correctBoundary = 1, stepSize = 1, startValue = 0, noiseSD = 0, decayAsymptote = 0.2, decayBeta = 0.0, loops = 200) {

  df.dist <- pCrossBoundaryDist(overlap, boundary, stepSize = stepSize, startValue = startValue, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops = loops)

  tmpMomentsA <- getMomentsOfDRM(df.dist$pCrossA, df.dist$samples)
  tmpMomentsA$overlap <- overlap
  tmpMomentsA$boundary <- boundary
  tmpMomentsA$startValue <- startValue
  tmpMomentsA$noiseSD <- noiseSD
  tmpMomentsA$decayAsymptote <- decayAsymptote
  tmpMomentsA$decayBeta <- decayBeta
  tmpMomentsA$correct <- ifelse(correctBoundary == 1, TRUE, FALSE)

  tmpMomentsB <- getMomentsOfDRM(df.dist$pCrossB, df.dist$samples)
  tmpMomentsB$overlap <- overlap
  tmpMomentsB$boundary <- -1*boundary
  tmpMomentsB$startValue <- startValue
  tmpMomentsB$noiseSD <- noiseSD
  tmpMomentsB$decayAsymptote <- decayAsymptote
  tmpMomentsB$decayBeta <- decayBeta
  tmpMomentsB$correct <- ifelse(correctBoundary == -1, TRUE, FALSE)

  tmpOut <- rbind(tmpMomentsA, tmpMomentsB)
  return(tmpOut)
}
