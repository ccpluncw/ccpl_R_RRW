#' This function takes the raw data and the rrwModelList as input, collapses the data, codes the columns, gets the parameter information, runs the smartGridSearch, and outputs the smartGridSearch results.
#'
#' Function that the raw data and the rrwModelList as input, collapses the data, codes the columns, gets the parameter information, runs the smartGridSearch, and outputs the smartGridSearch results.
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
#' @param numLoops An integer that specifies the number of loops with different parameter values to be searched in each run of the smartGridSearch before the best fits for that run are identified. DEFAULT = 200.
#' @param numIntervals An integer that specifies the number of intervals that the range between the upper and lower bounds of each parameter space will be segmented into.   DEFAULT = 25.
#' @param optParamListN An integer that specifies the number of best fit runs that will be used to extract the new bounds for next set of grid search loops.  So, if optParamListN == 10, then the parameter values from the 10 best runs will be used to refine the upper and lower bounds of for each parameter. DEFAULT = 10.
#' @param optBoundLoops An integer that specifies the number of times to re-simulate the best fit run in order to get a mean and SD for the minimized fit statistic. If the mean fit statistic is no longer the best fit, the program will move to the next best fit run, and so on until the it identifies the best mean fit run. DEFAULT = 10.
#' @param multicore An boolean that specifies whether the program should be run using multiple processing cores. DEFAULT = FALSE.
#' @param multicorePackages An vector of strings specifying the package names that are used in "fn." This needs to be used when multicore = TRUE. DEFAULT = NULL.
#' @param fileTag A string that is appended to the name of files to identify the analysis and experiment. The default is NULL, whereby the filetag will just be based on a timestamp.
#' @param numSimsToAverage A number specifying how many simulation runs are in the dataset and should be averaged together to get the final fit. Default is 10
#''
#' @return Several files are saved, including <fileTag>GridOut.txt (which contains the best 10 parameter fits), <fileTag>rrwData.txt (which contains the inputData), <fileTag>FitResults.txt (which contains the final parameter fits), <fileTag>Fitted.txt (which contains the input data + simulation fit data). The function returns a list containing (1) a list of the parameter values and fit statistics [runStats], and (2) a dataframe containing the input data with the best fit prediction based on the parameter values [df.fitted].
#' @keywords RRW smartGridSearch run
#' @export
#' @examples rrwRunSmartGridSearch(tmp.df, tmpModelList, minN = 40, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct01", correctVals = c(1,0), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = TRUE, numLoops = 2500, numIntervals = 100, optParamListN = 10, optBoundLoops = 10, multicore = TRUE, multicorePackages = c('RRW'), numSimsToAverage = 40, fileTag = NULL)

rrwRunSimulatedAnnealing <- function(data, rrwModelList, minN = 20, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct",correctVals = c(TRUE, FALSE), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = FALSE, pars.n = NULL, fileTag = NULL, numSimsToAverage = 10, finalNumSimsToAverage = 20, optim_sa_contol = list(t0 = 2000, nlimit = 200, r = 0.5), hoursToUpdate = 1, ...) {

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
  df.raw <- filterRRWdataMinN(df.raw, grpCols, "rt", "pHit", "n",dataOverlapCol,  minN)

  #output the rrwData to a file
  filename <- paste(fileTag,"rrwData.txt")
  write.table(df.raw, filename, col.names=T, row.names=F, quote=F, sep="\t")

  #get parameter information
  pList <- rrwGetModelParameterSpecs(rrwModelList)

  #add the necessary parameters to the "otherParameterList"
  otherParamList <- list(dataOverlapCol = dataOverlapCol, RwSamplesCol = RwSamplesCol, dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, numSimsToAverage = numSimsToAverage)

  #move the bounds, columns, and pars.n from pList to the appropriate format
  pIn <- pList$pU
  pU <- pList$pU
  pL <- pList$pL
  pUpper <- vector()
  pLower <- vector()
  pStart <- vector()
  df.parInfo <- NULL
  parClassNames <- names(pU)
  idx.1 <- 1

  for(pcn in parClassNames) {
    tmpParName <- names(pU[[pcn]])
    for(ppn in tmpParName) {
      if(pU[[pcn]][[ppn]] != pL[[pcn]][[ppn]]) {
        pUpper[idx.1] <- pU[[pcn]][[ppn]]
        pLower[idx.1] <- pL[[pcn]][[ppn]]
        pStart[idx.1] <- mean(c(pU[[pcn]][[ppn]], pL[[pcn]][[ppn]]))
 			  df.parInfo.tmp <- data.frame(parClass = pcn, parName = ppn)
 			  df.parInfo <- rbind(df.parInfo, df.parInfo.tmp)
 			  pIn[[pcn]][[ppn]] <- NA
        idx.1 <- idx.1 + 1
      }
    }
  }

  if(is.null(pars.n)) {
    otherParamList <- c(otherParamList, pList$pars.n)
  } else {
    otherParamList <- c(otherParamList, pars.n = pars.n)
  }

  otherParamList <- c(otherParamList, pList$pCols)
  otherParamList <- c(otherParamList, list(data = df.raw))

  rrwSim <- function(params_to_optimize) {

    parList <- list()

    numPars <- nrow(df.parInfo)

		for(i in 1:numPars) {
			pcn <- df.parInfo$parClass[i]
			ppn <- df.parInfo$parName[i]
		 	pIn[[pcn]][[ppn]] <- params_to_optimize[i]
		}

    if(!is.null(otherParamList)) {
      parList <- c(pIn, otherParamList)
    }

    outFit <- do.call(assessAveRRWfit,parList)

    ### Update every hour
    current_time <- Sys.time()
    time_since_last_update <- difftime(current_time, last_update_time, units = "hours")
    if (time_since_last_update >= hoursToUpdate) {
      cat("Hourly Update at:", format(current_time), "- Fit Statistic:", outFit, "\n")
      last_update_time <<- current_time
    }

    return(outFit)
  }

  rrwEnd <- function(params_to_optimize) {
    parList <- list()

    numPars <- nrow(df.parInfo)

    for(i in 1:numPars) {
      pcn <- df.parInfo$parClass[i]
      ppn <- df.parInfo$parName[i]
      pIn[[pcn]][[ppn]] <- params_to_optimize[i]
    }

    if(!is.null(otherParamList)) {
      parList <- c(pIn, otherParamList)
    }

    #create an input list to run "getMeanRRWfit" tht will run the simulation with the best fit
    #parameters multiple times and calculate the average of the runs.  The average is important to
    parList$minimizeStat <- NULL
    #add the output filenames
    parList$sinkFilename <- paste(fileTag,"FitResults.txt")
    #add the other paramters
    #How many simulations do you run and average together?
    parList$numSimsToAverage = finalNumSimsToAverage

    #save the Parameter values as an r data format
    filename <- paste(fileTag,"x.in.RData")
    save(parList, file=filename)
    #run the RRW to get the fit given the parameters output by the smartGridSearch
    run.list <- do.call(getMeanRRWfit, parList)
    df.fitted <-run.list$df.fitted
    #save the fit data
    filename <- paste(fileTag,"Fitted.txt")
    write.table(df.fitted, filename, col.names=T, row.names=F, quote=F, sep="\t")

    #add the minimization statistic to the output file
    sink(parList$sinkFilename, append = TRUE)
      cat("\n Minimization Statistic = ", minimizeStat, "\n")
    sink(NULL)

    return(run.list)
  }

  #get the start time
  start <- Sys.time()
  last_update_time <- start

  sa_result <- optimization::optim_sa(
    start = pStart,
    fun = rrwSim,
    lower = pLower,
    upper = pUpper,
    control = optim_sa_contol
  )

  ## get the end time
  end <- Sys.time()
  # see how long the analysis took to run
  diff <- end-start
  print(diff)

  outList <- rrwEnd(sa_result$par)

  outList[["sa_result"]] <- sa_result

  return(outList)
}
