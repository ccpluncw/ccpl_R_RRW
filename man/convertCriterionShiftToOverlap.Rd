% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertCriterionShiftToOverlap.r
\name{convertCriterionShiftToOverlap}
\alias{convertCriterionShiftToOverlap}
\title{Function that converts a evaluation criterion shift (in SD units) to an overlap. This is essentially identical to "convertValueShiftToOverlap()"}
\usage{
convertCriterionShiftToOverlap(
  currentOverlap,
  criterionShiftInSDunits,
  allowExactZeroToShift = T
)
}
\arguments{
\item{currentOverlap}{A number that is the current overlap of the two distributions before the mean shift.}

\item{criterionShiftInSDunits}{A number that is the evaluation criterion shift shift SD units.}

\item{allowExactZeroToShift}{A boolean that specifies whether to allow an exact zero overlap to shift. You may want to set this to FALSE if you beleive that an overlap of Zero indicated distributions that are so far appart that any mean shift if unlikely to change that overlap. When it is set to TRUE, the overlap of Zero is set to 0.0001 (allowing a shift to occur). This only applies to "exact zeros." When this function reduces an overlap to Zero, that reduction size is known and will not require a correction. DEFAULT = TRUE.}
}
\value{
A number that represents the new overlap after the mean shift is taken into account
}
\description{
This function converts a evaluation criterion shift (in SD units) to an overlap. To do so, it uses an equation fit to normal distribution changes - which will fit the average distibution shape (if not all individual distributions)
}
\examples{
convertCriterionShiftToOverlap (0.5, 1)
}
\keyword{change}
\keyword{mean}
\keyword{overlap}
\keyword{shift}
\keyword{value}
\keyword{vc}
