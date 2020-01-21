#' This function simulates the RRW for a single set of parameters across a single overlap.
#'
#' This function simulates the RRW for a single set of parameters across a single overlap. It returns a dataframe that describes the probability of crossing each boundary for each sample size.  It presents sample sizes from 1 to max, where max is when pCrossA + pCrossB = 1.0.  If pCrossA + pCrossB < 1.0 and pCrossA + pCrossB does not change for the last 10 percent of the samples, then boundary = boundary - 0.1.  Thus, the boundaries move closer together slowly if the walk is not crossing a boundary.
#' @param overlap A number between 0-1 representing a single distributional overlap.
#' @param boundary A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param stepSize A number that specifies the granularity that the simulation will generate the probability crossing a boundary.  stepSize = 1 is the Default and most granular.
#' @param startValue A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  startValue = 0 is the default and represents an unbiased start point.
#' @param noiseSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. noiseSD = 0 is the default and represents no noise being added to each step.
#' @param decayAsymptote A positive number representing the asymptote value in the Information Accrual Bias (IAB). decayAsymptote = 0.2 is the default.
#' @param decayBeta A signed number representing the beta value in the Information Accrual Bias (IAB). decayBeta = 0 is the default and represents no IAB.
#' @param loops A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @return A dataframe that contains the boundary; samples; overlap; the p(D1 > D2) derived from overlap; pCrossA the probability of the RRW with this number of samples crossing the possitive boundary; pCrossB the probability of the RRW with this number of samples crossing the negative boundary.
#' @keywords RRW random walk simulation pCross Boundary
#' @export
#' @examples pCrossBoundaryDist (overlap=0.1, boundary=14, startValue=0.1, loops = 400)

pCrossBoundaryDist <- function (overlap, boundary, stepSize = 1, startValue = 0, noiseSD = 0, decayAsymptote = 0.2, decayBeta = 0.0, loops = 200) {

  df.dist <- NULL
  i <- 1
  step <- 1
  done <- FALSE
  while(done == FALSE) {
    df.tmp <- RRW(overlap, boundary, step, startValue = startValue, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops=loops)
    df.dist <- chutils::ch.rbind(df.dist, df.tmp)

    if(i > 10) {
      done <- ifelse( df.dist$pCrossA[i] + df.dist$pCrossB[i] >= 0.999, TRUE, FALSE)
    }

    if(i > 100) {
      #if there has been little change in the last 10% of trials
      if ( df.dist$pCrossA[length(df.dist$pCrossA)] <= mean(df.dist$pCrossA[(length(df.dist$pCrossA) - length(df.dist$pCrossA)/10):length(df.dist$pCrossA)]) ) {
        #and boundary is greater than 0
        if(boundary > 0) {
          #then reduce the boundary
          boundary <- boundary - 0.1
        }
      }
    }

    i <- i + 1
    step <- step + stepSize
  }

  df.dist$sumAB <- df.dist$pCrossA + df.dist$pCrossB

  return(df.dist)

}
