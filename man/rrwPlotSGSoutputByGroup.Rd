% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rrwPlotSGSoutputByGroup.r
\name{rrwPlotSGSoutputByGroup}
\alias{rrwPlotSGSoutputByGroup}
\title{This function plots the rrw data fits from the output of getMeanRRWfit() and the rrwModelList.}
\usage{
rrwPlotSGSoutputByGroup(
  df.RRWout,
  grpCols = NULL,
  dataRtCol = "rt",
  dataPhitCol = "pHit",
  rtFitCol = "rtFit",
  pHitFitCol = "pCross",
  correctCol = "correct",
  overlapCol = "overlap",
  fileTag = NULL,
  numSimsToPlot = 10
)
}
\arguments{
\item{df.RRWout}{This is a dataframe containing the output from getMeanRRWfit() or rrwRunWithFixedParameters () or rrwRunSmartGridSearch (). It contains the rrwData and the simulation fits.}

\item{grpCols}{A vector of strings containing the grouping columns. If NULL this argument will be ignored.  DEFAULT = NULL}

\item{dataRtCol}{A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"}

\item{dataPhitCol}{A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"}

\item{rtFitCol}{A string that identifies the name of the column in data that contains the RRW fit of the RTs for the specific overlap/correct/condition combination. Default is "rt"}

\item{pHitFitCol}{A string that identifies the name of the column in data that contains the RRW fit of the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pCross"}

\item{correctCol}{a boolean that specifies the name of the column that specifies if response was correct or incorrect.  It must be coded as a boolean, whereby correct == TRUE, and incorrect == FALSE. The default is "correct"}

\item{overlapCol}{A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"}

\item{fileTag}{A string that is appended to the name of files to identify the analysis and experiment. The default is NULL, whereby the filetag will just be based on a timestamp.}

\item{numSimsToPlot}{A number specifying how many simulation runs are in the dataset and should be plotted. Default is 10.
'}
}
\value{
Several files containing plots are saved.  The plots are duplicated: one version with each simulation plus the average and one with just the average.
}
\description{
Function that plots the rrw data fits from the output of getMeanRRWfit() and the rrwModelList.
}
\examples{
rrwPlotSGSoutput(df.fitted, "rt", "pHit", "rtFit", "pCross", "correct", "overlap", fileTag = NULL, numSimsToPlot = 40)
}
\keyword{RRW}
\keyword{plots}
