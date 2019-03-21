#' This function simulates the RRW to build the response distribution.
#'
#' This function simulates the RRW to build the response distribution.  It does so incrementally, one numSamples at a time.  Specifically, for a specific overlap and all the usual parameters of the RRW model, it will assess the probability of crossing the boundary given a specific number of samples ("samples"). The function pCrossBoundaryDist() calls RRW multiple times to build the entire distribution for a specific overlap and all the usual parameters.
#'
#' @param overlap A number between 0-1 representing a single distributional overlap.
#' @param boundary A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param samples The number of samples in this particular run.  Generally, this will be input from 1 up in steps of 1 until all of the samples cross a boundary.
#' @param startValue A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  startValue = 0 is the default and represents an unbiased start point.
#' @param noiseSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. noiseSD = 0 is the default and represents no noise being added to each step.
#' @param decayAsymptote A positive number representing the asymptote value in the Information Accrual Bias (IAB). decayAsymptote = 0.2 is the default.
#' @param decayBeta A signed number representing the beta value in the Information Accrual Bias (IAB). decayBeta = 0 is the default and represents no IAB.
#' @param loops A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#'
#' @return A dataframe that contains the boundary; samples; overlap; the p(D1 > D2) derived from overlap; pCrossA the probability of the RRW with this number of samples crossing the possitive boundary; pCrossB the probability of the RRW with this number of samples crossing the negative boundary.
#' @keywords RRW random walk simulation moments overlap
#' @export
#' @examples RRW (overlap = 0.1, boundary=14, samples = 10, startValue=0.1, loops = 400)
#' @details The Robust Random Walk (RRW) sequential sampling procedure is different from traditional sequential sampling procedures in two ways.  First, the RRW does not estimate drift rate.   Rather, it takes Overlap as the input. Overlap equals the probability that a sampled value from Distribution 1 (D1) would be greater than a sampled value from Distribution 2 (D2).  We convert this probability into a measure of overlap, such that,
#'
#'              Overlap=1-(abs(p(D1 > D2)-0.5)/0.5)
#'
#' Here, Overlap = 0 indicates that p(D1 > D2) = 0 or 1, and Overlap = 1 indicates that p(D1 > D2) = 0.5.  Furthermore, when Overlap<1, we note which item has the higher p(D1 > D2).
#'
#' 	The second way that RRW differs from traditional sequential sampling procedures is that it is non-parametric. The most common sequential sampling procedures assume that the latent psychological distributions driving decision choice have a specific distributional shape, such as Gaussian or Laplace.  The RRW, in contrast, makes no such assumptions.  Rather than estimate Overlap as a free parameter, the RRW takes Overlap as a predictor variable.  Because Overlap is equivalent to the p(D1 > D2), the RRW uses p(D1 > D2) to drive the accumulation of evidence.  Specifically, each step, i, of the RRW equals:
#'
#'            Ev = ((1 if p(D1 > D2) or -1 if p(D1 < D2)) + e) * IAB
#'
#'  Where,
#'
#'            e ~ N(0,noiseSD)
#'
#'  And
#'
#'            decay= (1- decayAsymptote )*exp⁡(decayBeta *i)+ decayAsymptote
#'
#'  Briefly, with each sample of the random walk, EVi is assigned a 1 with the probability of p(D1 > D2) and a -1 with the probability p(D1 < D2).  These probabilities are taken directly from the magnitude estimation data.  To that, we add noise sampled from a normal distribution with a mean of 0 and a standard deviation of noiseSD. If one wished to add noise, noiseSD will be greater than 0 and will be a free parameter in the model.  IAB stands for "information accrual bias."  The IAB is an exponential function weights the influence of the sample, EV, by the time since the start of the process.  If decayBeta is negative, then the early samples have more influence on the random walk than the later samples.  This corresponds to a primacy of information bias. Conversely, if decayBeta is positive, then the later samples have more influence on the random walk than the early samples.  This corresponds to a recency of information bias.  The decayAsymptote and decayBeta work in unison to describe the shape and strength of this bias.  If either decayBeta = 0 or decayAsymptote = 1, then IAB = 1, so there is not information accrual bias.  In general, as decayBeta moves away from 0 (positive or negative), the influence of IAB increases. When decayBeta < 0, dA acts as a boundary, such that IAB will not decrease below its value.  When decayBeta > 0, decayAsymptote acts as a moderator, such that as decayAsymptote approaches 1 it will reduce the influence of IAB.  Both dA and dB can be free parameters.  We often fix decayAsymptote = 0.2, and make dB a free parameter.
#'
#'	In addition to the free parameters, the we also include the following free parameters:
#'
#'	boundary distance from a 0 startValue, b. Within the RRW, the distance between boundaries is symmetric around a startValue of 0. Within the RRW, the positive boundary corresponds to the correct response to the higher valued option, whereas the negative boundary corresponds to the incorrect response to the higher valued option.
#'
#'	The RRW includes the startValue as a signed proportion of the boundary distance from 0, startValue.  If startValue > 0, it represents a bias towards the positive boundary, with the starting position equaling startValue*BOUNDARY. If startValue < 0, it represents a similar bias towards the negative boundary.
#'
#'	Finally, we include the non-decision time, Ter, to represent what is generally measured as the intercept of the function. It is hypothesized to include all ancillary processes not associated with the variable of interest (e.g., Psychological Value).
#'
#'  The RRW also has a scaling parameter that scales the number of steps into RT.  We do not count this scaling parameter as a free parameter because it is a simple linear transformation of one relatively arbitrary scale into another.

#'

RRW <- function (overlap, boundary, samples, startValue = 0, noiseSD = 0, decayAsymptote = 0.2, decayBeta = 0.0, loops = 200) {

  ## convert overlap into a p value
  p <- overlapToP(overlap)

  if(abs(startValue) >= 1) {
    print("startValue = ")
    print(startValue)
    stop("startValue must be a value between -1 and 1. startValue * boundary = the starting point of the random walk ")
  } else {
    startPoint <- startValue * boundary
  }

  ## set the stepSize (vec) and probabilities for the Non-parametric walk (prob)
  vec <- c(1,-1)
  prob <- c(p, (1-p))

  ### run a loop to simulate the random walk
  sumCrossA <- 0
  sumCrossB <- 0
  stepsSeq <- c(1:samples)
  for(i in 1:loops) {

    ### this is the main call
    ### collect samples, and add noise if noiseSD is > 0
    outVec <-  sample(vec, samples, prob=prob, replace = TRUE) +  rnorm(samples, 0, abs(noiseSD))

    ### add decay over time.  When decayBeta == 0, then no decay; decayBeta < 0, less info over time; decayBeta > 0, more info over time
    decayC <- (1-decayAsymptote)*exp(decayBeta*stepsSeq)+decayAsymptote
    outVec <- outVec * decayC

    ### add the startValue to the first run
    outVec <- append(outVec, startPoint, 0)

    ## identify if/when the walk crosses a boundary
    outVecCS <- cumsum(outVec)
    maxI <- ifelse(max(outVecCS) >= boundary, min(which(outVecCS >= boundary)), -1)
    minI <- ifelse(min(outVecCS) <= (-1*boundary), min(which(outVecCS <= (-1*boundary))), -1)

    ## identify which boundary was crossed
    if(maxI > 0) {
      if(minI > 0) {
        if(maxI < minI) {
          sumCrossA <- sumCrossA + 1
        } else {
          sumCrossB <- sumCrossB + 1
        }
      } else {
        sumCrossA <- sumCrossA + 1
      }
    } else {
      if(minI > 0) {
        sumCrossB <- sumCrossB + 1
      }
    }
  }

  #output data
  return(data.frame(boundary = boundary, samples = samples, overlap = overlap, p = p, pCrossA = sumCrossA/loops, pCrossB = sumCrossB/loops))
}
