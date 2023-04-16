#' This function runs a set of fixed parameters (specified in a rrwModelList) to generate the fits to the data. It takes the raw data and the rrwModelList as input, collapses the data, codes the columns, gets the parameter information, runs the getMeanRRWfit, and outputs the smartGridSearch results. rrwModelList should only have fixed parameter values (upperBound==lowerBound).
#'
#' Function that runs a set of fixed parameters (specified in a rrwModelList) to generate the fits to the data.
#' @param data This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param minN A a numeric value that specifies the minimum number of trials necessary for a condition to be included in the analysis.  If the condition has less than minN trials, then that condition is assigned and NA.  Default=20
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param RwSamplesCol A string that identifies the name of the column in the RRW simulations that contains the summary of the samples that you want to use as a simulation for RT.  The possible columns are: "Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile).  The default is "Q50" because the median is more robust than the mean.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param correctCol a string that specifies the name of the column that specifies if the participant chose the item with the greatest value distribution (correct) or if they did not (incorrect).
#' @param correctVals a vector of two values that specifies the "correct" value (index 1) and the "incorrect" value (index 2). e.g, c("yes", "no"). DEFAULT = c(TRUE, FALSE)
#' @param loopsPerRWstep A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @param minimizeStat A string that specifies which statistic to minimize when optimizing the model fit.  The options are: "BIC" , "AIC" , or "R_Square". Default is "BIC".
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.
#' @param pars.n The number of free parameters. When NULL, the program will attempt to calculate the number of free parameters from the input. Default = NULL.
#' @param fileTag A string that is appended to the name of files to identify the analysis and experiment. The default is NULL, whereby the filetag will just be based on a timestamp.
#' @param numSimsToAverage A number specifying how many simulation runs are in the dataset and should be averaged together to get the final fit. Default is 10.
#''
#' @return Several files are saved, including <fileTag>rrwData.txt (which contains the inputData), <fileTag>FitResults.txt (which contains the final parameter fits), <fileTag>Fitted.txt (which contains the input data + simulation fit data). The function returns a list containing (1) a list of the parameter values and fit statistics [runStats], and (2) a dataframe containing the input data with the best fit prediction based on the parameter values [df.fitted].
#' @keywords RRW fixed parameter run getMeanRRWfit
#' @export
#' @examples rrwRunWithFixedParameters(tmp.df, tmpModelList, minN = 40, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct01", correctVals = c(1,0), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = TRUE, numSimsToAverage = 40, fileTag = NULL)

rrwRunWithFixedParameters <- function(data, rrwModelList, minN = 20, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct",correctVals = c(TRUE, FALSE), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = FALSE, pars.n = NULL, fileTag = NULL, numSimsToAverage = 10, ...) {

  #if fileTag does not exist, create one based on the timestamp
  if(is.null(fileTag)) {
    fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), "RRW", sep="_")
  }

  #get grouping columns from modelList
  grpCols <- getGroupColsFromModelList(rrwModelList)
  #collapse data
  df.raw <- rrwRawDataToRRWFormatData(data, grpCols, dataRtCol, correctCol,correctVals, dataOverlapCol)
  #add coded columns
  df.raw <- rrwGetAllCodeColumns(df.raw, rrwModelList)

  #remove low n conditions
  #df.raw$rt <- ifelse(df.raw$n < minN, NA, df.raw$rt)
  df.raw <- filterRRWdataMinN(df.raw, grpCols, "rt", "pHit", "n",dataOverlapCol,  minN)
  #output the rrwData to a file
  filename <- paste(fileTag,"rrwData.txt")
  write.table(df.raw, filename, col.names=T, row.names=F, quote=F, sep="\t")

  #get parameter information
  pList <- rrwGetModelParameterSpecs(rrwModelList)

  #add the necessary parameters to the "otherParameterList"
  otherParamList <- list(dataOverlapCol = dataOverlapCol, RwSamplesCol = RwSamplesCol, dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit)

  #move the bounds, columns, and pars.n from pList to the appropriate format
  pU <- pList$pU
  pL <- pList$pL
  pI <- pList$pI
  otherParamList <- c(otherParamList, pList$pCols)
  if(is.null(pars.n)) {
    otherParamList <- c(otherParamList, pList$pars.n)
  } else {
    otherParamList <- c(otherParamList, pars.n = pars.n)
  }
  otherParamList <- c(otherParamList, list(data = df.raw))

  #create an input list to run "getMeanRRWfit" tht will run the simulation with the best fit
  #parameters multiple times and calculate the average of the runs.  The average is important to
  #smooth out the random fluctuations of the stochastic process.

  #Since we are using a modelList,  upperBound == lowerBound when a parameter is fixed, so we will input the upperBound here
  x.in <- as.list(pU)
  #NULL out the unused parameters
  otherParamList$minimizeStat <- NULL
  #add the output filenames
  x.in$sinkFilename <- paste(fileTag,"FitResults.txt")
  #add the other paramters
  x.in <- c(x.in, otherParamList)
  #How many simulations do you run and average together?
  x.in$numSimsToAverage = numSimsToAverage

  #save the Parameter values as an r data format
  filename <- paste(fileTag,"x.in.RData")
  save(x.in, file=filename)
  #run the RRW to get the fit given the parameters output by the smartGridSearch
  run.list <- do.call(getMeanRRWfit, x.in)
  df.fitted <-run.list$df.fitted
  #save the fit data
  filename <- paste(fileTag,"Fitted.txt")
  write.table(df.fitted, filename, col.names=T, row.names=F, quote=F, sep="\t")

  #add the minimization statistic to the output file
  sink(x.in$sinkFilename, append = TRUE)
    cat("\n Minimization Statistic = ", minimizeStat, "\n")
  sink(NULL)

  return(run.list)
}
