% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rrwRunStatsToDataframe.r
\name{rrwRunStatsToDataframe}
\alias{rrwRunStatsToDataframe}
\title{This function converts the runStats list (which is the output of the rrw fit functions) into a single row of a dataframe.}
\usage{
rrwRunStatsToDataframe(runStats)
}
\arguments{
\item{runStats}{A list that contains the output of a single run of the kFold cross validation of the RRW.  It is output by rrw fit functions and used in the rrwKfoldStatsToDataframe, etc.
'}
}
\value{
a single row of a dataframe that contains the parameter values and fit statistics in each column.
}
\description{
Function that converts the runStats list into a single row of a dataframe.
}
\examples{
rrwKfoldStatsToDataframe (runStats)
}
\keyword{RRW}
\keyword{kFold}
\keyword{rrw}
\keyword{runStats}
