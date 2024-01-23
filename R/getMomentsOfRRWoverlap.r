#' This function simulates the RRW for a single set of parameters across all values of overlap.
#'
#' This function simulates the RRW for a single set of parameters across all values of overlap.
#' @param overlapSeq A vector of numbers between 0-1 representing distributional overlaps.
#' @param boundary A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param boundaryChangeSensitivity A positive number that adjusts the likelihood that the boundary will be reduced because there it is taking too many samples before reaching a threshold. Essentially, this models a decision-maker reducing the amount of information that they need before a response because it is taking to long.  The boundary reduces when there is little change in the likelihood of a response after N number of steps.  The routine makes this determination as a function of the size of the unaltered boundary (boundaryChangeSensitivity * boundary). Small values (e.g., 0.1), result in an impatient person (relatively early chnage in the boundary values).  Large values (e.g., 0.9), result in a patient person (a boundary that is relatively resistent to change.  A value of 0 is an exception: it will result in no boundary reduction under any circumstances. Impatience increases the likelihood of an error, so it increases the average time of error responses.  Default is 0.25.
#' @param Ter A positive number that adjusts the time of encoding/response (non-decision time). When Ter = 0, a single, best fitting Ter will be calculated. Use this parameter when you have several conditions and you hypothesize that some conditions should have a larger Ter than others.  For example, this may model longer encoding process when different numbers of items are on the screen for different conditions. The hypothetical fastest condition should have Ter = 0, which simply indicates that the RRW program will find the best fit Ter. The slower conditions will have larger Ter values and will be fit ralative to the fastest condition. The change in samples resulting from the Ter effect is a function of the size of the boundary (Ter * boundary). Default is 0.0. This is still experimental.
#' @param stepSize A number that specifies the granularity that the simulation will generate the probability crossing a boundary.  stepSize = 1 is the Default and most granular.
#' @param startValue A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  startValue = 0 is the default and represents an unbiased start point.
#' @param noiseSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. noiseSD = 0 is the default and represents no noise being added to each step.
#' @param decayAsymptote A positive number representing the asymptote value in the Information Accrual Bias (IAB). decayAsymptote = 0.2 is the default.
#' @param decayBeta A signed number representing the beta value in the Information Accrual Bias (IAB). decayBeta = 0 is the default and represents no IAB.
#' @param loops A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @param correctBoundary Either 1 or -1.  If correctBoundary = 1, then the positive boundary is the boundary for a correct response.  If correctBoundary = -1, then the negative boundary is the boundary for a correct response.
#' @param progress TRUE or FALSE that specifies whether to present a progress bar.  Default is FALSE.
#' @return A dataframe that contains the moments of the simmulated values from the model ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the probability of crossing each boundary); correct (TRUE or FALSE)), along with the parameters of the model.
#' @keywords RRW random walk simulation moments overlap
#' @export
#' @examples getMomentsOfRRWoverlap (seq(0,1,0.1), boundary=14, startValue=0.1, loops = 400)

getMomentsOfRRWoverlap <- function (overlapSeq, boundary, boundaryChangeSensitivity = 0.25, Ter = 0, stepSize = 1, startValue = 0, noiseSD = 0, decayAsymptote = 0.2, decayBeta = 0.0, loops = 200, correctBoundary = 1, progress = FALSE) {

  if(progress) {
    pb <- txtProgressBar(1, length(overlapSeq), 1, style = 3, width=length(overlapSeq))
  }
  df.rwAll <- NULL
  for(i in 1: length(overlapSeq)) {
      df.rwTmp <- getMomentsOfRRW(overlap = overlapSeq[i], boundary, boundaryChangeSensitivity = boundaryChangeSensitivity, Ter = Ter, stepSize = stepSize, startValue = startValue, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops=loops, correctBoundary = correctBoundary)

      df.rwAll <- chutils::ch.rbind(df.rwAll,df.rwTmp)
      if(progress) {
        setTxtProgressBar(pb, i)
      }
  }

  return(df.rwAll)

}
