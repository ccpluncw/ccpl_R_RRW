% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rrwCodeColumn.r
\name{rrwCodeColumn}
\alias{rrwCodeColumn}
\title{This  function uses the rrwModelList to add the effect and dummy variables to the RRW data (one column at a time)}
\usage{
rrwCodeColumn(data, rrwModelList)
}
\arguments{
\item{data}{This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.}

\item{rrwModelList}{A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
'}
}
\value{
the data is returned with the added dummy/effect coded column
}
\description{
Function that uuses the rrwModelList to add the effect and dummy variables to the RRW data (one column at a time)
}
\examples{
rrwCodeColumn (data, modelList)
}
\keyword{RRW}
\keyword{code}
\keyword{column}
\keyword{data}
\keyword{dummy}
\keyword{effect}
\keyword{rrw}
