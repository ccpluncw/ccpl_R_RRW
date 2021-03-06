% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRowParameterValue.r
\name{getRowParameterValue}
\alias{getRowParameterValue}
\title{This function outputs the parameter value for a specific row of the data.}
\usage{
getRowParameterValue(data.row, parCol = NULL, parVec)
}
\arguments{
\item{data.row}{This is a single row of the  dataframe that  contains (at least) the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.}

\item{parCol}{This is a vector of strings that identify the names of the columns in data that contain the effect codes for the specific parameter entered here. Default  = NULL (e.g., no effect codes for this parameter)}

\item{parVec}{A vector of number(s) specifying the values for this specific paramter.}
}
\value{
The value(s) of the paramter for the condition specified in data.row.
}
\description{
This function outputs the parameter value for a specific row of the data. Because models are effect coded, each parameter value for the data row is often the sum of the effects.  This function calculates and outputs that sum.
}
\examples{
getRowParameterValue (data[1,], "boundaryCols", b)
}
\keyword{RRW}
\keyword{codes}
\keyword{effect}
\keyword{parameters}
\keyword{random}
\keyword{walk}
