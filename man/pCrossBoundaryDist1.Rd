% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pCrossBoundaryDistold1.r
\name{pCrossBoundaryDist1}
\alias{pCrossBoundaryDist1}
\title{This function simulates the RRW for a single set of parameters across a single overlap.}
\usage{
pCrossBoundaryDist1(
  overlap,
  boundary,
  stepSize = 1,
  boundaryChangeSensitivity = 0.25,
  startValue = 0,
  noiseSD = 0,
  decayAsymptote = 0.2,
  decayBeta = 0,
  loops = 200
)
}
\arguments{
\item{overlap}{A number between 0-1 representing a single distributional overlap.}

\item{boundary}{A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value}

\item{stepSize}{A number that specifies the granularity that the simulation will generate the probability crossing a boundary.  stepSize = 1 is the Default and most granular.}

\item{boundaryChangeSensitivity}{A number that adjusts the likelihood that the boundary will be reduced because there it is taking too many samples before reaching a threshold. Essentially, this models a decision-maker reducing the amount of information that they need before a response because it is taking to long.  The boundary reduces when there is little change in the likelihood of a response after N number of steps.  The routine makes this determination as a function of the size of the unaltered boundary (boundaryChangeSensitivity * boundary). Small values (e.g., 0.1), result in an impatient person (relatively early chnage in the boundary values).  Large values (e.g., 0.9), result in a patient person (a boundary that is relatively resistent to change.  A value of 0 is an exception: it will result in no boundary reduction under any circumstances. Impatience increases the likelihood of an error, so it increases the average time of error responses.  Default is 0.25.}

\item{startValue}{A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  startValue = 0 is the default and represents an unbiased start point.}

\item{noiseSD}{A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. noiseSD = 0 is the default and represents no noise being added to each step.}

\item{decayAsymptote}{A positive number representing the asymptote value in the Information Accrual Bias (IAB). decayAsymptote = 0.2 is the default.}

\item{decayBeta}{A signed number representing the beta value in the Information Accrual Bias (IAB). decayBeta = 0 is the default and represents no IAB.}

\item{loops}{A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.}
}
\value{
A dataframe that contains the boundary; samples; overlap; the p(D1 > D2) derived from overlap; pCrossA the probability of the RRW with this number of samples crossing the possitive boundary; pCrossB the probability of the RRW with this number of samples crossing the negative boundary.
}
\description{
This function simulates the RRW for a single set of parameters across a single overlap. It returns a dataframe that describes the probability of crossing each boundary for each sample size.  It presents sample sizes from 1 to max, where max is when pCrossA + pCrossB = 1.0.  If pCrossA + pCrossB < 1.0 and pCrossA + pCrossB does not change for the last 10 percent of the samples, then boundary = boundary - 0.1.  Thus, the boundaries move closer together slowly if the walk is not crossing a boundary.
}
\examples{
pCrossBoundaryDist (overlap=0.1, boundary=14, startValue=0.1, loops = 400)
}
\keyword{Boundary}
\keyword{RRW}
\keyword{pCross}
\keyword{random}
\keyword{simulation}
\keyword{walk}
