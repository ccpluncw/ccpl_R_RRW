% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rrwAssessComparisonDistribitionFit.r
\name{rrwAssessComparisonDistribitionFit}
\alias{rrwAssessComparisonDistribitionFit}
\title{A function to assess how well a specific comparison distribution fit with choice data}
\usage{
rrwAssessComparisonDistribitionFit(
  df.choiceDat,
  df.valueDat,
  rdN,
  rdMean,
  rdSD,
  refDistGenFn = rnorm,
  item1cols,
  item2cols,
  respChoiceCol,
  respChoiceVal,
  RTcol,
  chanceThreshold = 0,
  lowRTquantileThreshold = 0,
  highRTquantileThreshold = 0,
  minOverlapN = 0,
  pars.n = 3,
  equalizeRTandPhit = TRUE,
  minimizeStat = "BIC",
  roundThreshold = 0.1,
  roundDirection = ceiling,
  numOverlapBins = 10,
  overlapNumRuns = 1000,
  overlapBootstrapMulticore = FALSE,
  useTwoParameterModel = FALSE,
  onlyMinimizeOnPhvoFit = FALSE,
  splitByRefDist = TRUE,
  fitByRefDist = FALSE,
  splitByCorrect = TRUE,
  combFun = ch.maxAveComb,
  ...
)
}
\arguments{
\item{df.choiceDat}{A dataframe with the choice data from a typical PVT experiment.}

\item{df.valueDat}{A dataframe with the raw values data from a typical values experiment.}

\item{rdN}{An integer specifying the number of observations in the reference distribution.}

\item{rdMean}{An number specifying the mean of the reference distribution (Gaussian).}

\item{rdSD}{An number specifying the standard deviation of the reference distribution (Gaussian).}

\item{refDistGenFn}{A function to generate the observation for the reference distribution.  This function must take the following arguments: n, mean, sd.  DEFAULT = rnorm.}

\item{item1cols}{a vector of strings that specifies the names of the columns in "df.choiceDat" that contains the the probes in Item 1.}

\item{item2cols}{a vector of strings that specifies the names of the columns in "df.choiceDat" that contains the the probes in Item 2.}

\item{respChoiceCol}{a string that specifies the name of the column in "df.choiceDat" that contains the the participant's response to the prompt - yes take action or no take no action.}

\item{respChoiceVal}{a vector of two values that specifies the choose Item1 option ("yes" take action in many morals experimants, thus saving Item1) value (index 1) and the choose Item1 option ("no" take no action in many morals experimants, thus saving Item2) value (index 2).}

\item{RTcol}{a string that specifies the name of the column in "df.choiceDat" that contains the RT for each trial.}

\item{chanceThreshold}{A number specifying the quantile that all subjects whose p(HVO) falls below it will be removed from the dataset. DEFAULT = 0 (none removed)}

\item{lowRTquantileThreshold}{A number specifying the quantile that, for each overlapRound, individual RTs that fall below it will be removed from the dataset. DEFAULT = 0.0 (none removed)}

\item{highRTquantileThreshold}{A number specifying the quantile that, for each overlapRound, individual RTs that fall above  it will be removed from the dataset. DEFAULT = 1.0 (none removed)}

\item{minOverlapN}{A number specifying the minimum number of responses an overlapRound condition must have to remain in the dataset. DEFAULT = 0 (none removed)}

\item{pars.n}{The number of free parameters. When NULL, the program will attempt to calculate the number of free parameters from the input. Default = 3.}

\item{equalizeRTandPhit}{A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.}

\item{minimizeStat}{A string that specifies which statistic to minimize when optimizing the model fit.  The options are: "BIC" , "AIC" , or "R_Square". Default is "BIC".}

\item{roundThreshold}{An integer that specifies the nearest interval that the overlaps should be rounded to. DEFAULT = 0.1 (round to the nearest 0.1)}

\item{roundDirection}{An option that specifies the rounding direction: ceiling (always round up), floor (always round down), or round (round to the nearest value, up or down). DEFAULT = ceiling}

\item{numOverlapBins}{An integer that specifies the number of bins the overlaps should be binned into. This is only used if roundThreshold = NULL.  DEFAULT = 10 (round overlaps so that there are 10 bins)}

\item{overlapNumRuns}{the number of runs to do in the bootstrap that calculates the overlap of the reference distribution and the choice itms. DEFAULT = 1000.}

\item{overlapBootstrapMulticore}{A boolean that specifies whether to run the bootstrap that calculates the value overlaps in multicore mode (in parallel).  This should be set to FALSE if you are using this function within the smartGridSearch. DEFAULT = FALSE.}

\item{useTwoParameterModel}{A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = FALSE.}

\item{onlyMinimizeOnPhvoFit}{A boolean that specifies whether to minimize the function only on the fit of the p(HVO) data (if set to TRUE).  If set to false, the function will be minimized on the basis of both the fit of the p(HVO) and RT data.  DEFAULT = FALSE.}

\item{splitByRefDist}{A boolean that specifies whether to group the dataset by those trial that are above the reference distribution and those that are below the reference distribution.  DEFAULT = TRUE.}

\item{fitByRefDist}{A boolean that specifies whether fit the p(HVO) and RT functions by the grouping created when splitByRefDist = TRUE. If splitByRefDist = FALSE, then this argument is irrelevant. Generally, you want splitByRefDist = TRUE and  fitByRefDist = FALSE to find a reference distribution that fits an unbiased observer for a particular condition in your experiment and then see how your manipulations influence either another reference distribution or the response biases in the RRW.  DEFAULT = FALSE.}

\item{splitByCorrect}{A boolean that specifies whether to use the "correct/incorrect" column as a grouping variable for the RT data when fitting the RT function. When set to FALSE, the RT function will be fit to the collapsed correct and incorrect responses. DEFAULT = TRUE.}

\item{combFun}{If there are multiple items contributing to a single distribution, this function describes how the values will be combined across items in both the X and Y distributions. The function must combine elements of a list that might be of different different lengths. The default just flattens the list into one large vector.  DEFAULT = ch.maxAveComb  (with probMax = 0.5)}
}
\value{
The minimization statistic for the fit of the model to the RT and pHit (proportion crossed each boundary) data.  This is the value that will be miniized by the optimization program.
}
\description{
This function assesses how well a specific comparison distribution fit with choice data. The comparison distribution is normal, with a mean and SD.
}
\examples{
rrwAssessComparisonDistribitionFit (df.choiceDat=moralsData, df.valueDat = allValues, rdN=100, rdMean=3, rdSD=1, item1cols = "Item1", item2cols = "Item2", respChoiceCol = "keyDef", respChoiceVal = c("Yes", "No"), RTcol = "res.rt",chanceThreshold = 0.5,  lowRTquantileThreshold = 0.025, highRTquantileThreshold = 0.975, minOverlapN = 20, pars.n = 3, equalizeRTandPhit = TRUE, minimizeStat = "BIC", roundThreshold = 0.1, roundDirection = ceiling, numOverlapBins = 10, overlapNumRuns = 1000, combFun = ch.maxAveComb, probMax = 0.5)
}
\keyword{assess}
\keyword{distribution}
\keyword{fit}
\keyword{reference}
