#' This function runs a k-fold cross validity on the RRW analysis.  Generally use this for the best fit model to determine how robust it is.
#'
#' Function that runs a k-fold cross validity on the RRW analysis.
#' @param kFold An integer that specifies the number of "folds" to cross validate the model with. DEFAULT = 10.
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
#' @return Several files are saved, including <fileTag>GridOut.txt (which contains the best 10 parameter fits), <fileTag>rrwData.txt (which contains the inputData), <fileTag>FitResults.txt (which contains the final parameter fits), <fileTag>Fitted.txt (which contains the input data + simulation fit data). The function returns a list containing (1) a dataframe of the parameter values and fit statistics for all the k-fold runs.
#' @keywords RRW smartGridSearch run
#' @export
#' @examples kfoldRRWRunSmartGridSearch(kfold = 10, tmp.df, tmpModelList, minN = 40, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct01", correctVals = c(1,0), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = TRUE, numLoops = 2500, numIntervals = 100, optParamListN = 10, optBoundLoops = 10, multicore = TRUE, multicorePackages = c('RRW'), numSimsToAverage = 40, fileTag = NULL)

kfoldRRWRunSmartGridSearch <- function(kFold = 10, data, rrwModelList, minN = 20, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct",correctVals = c(TRUE, FALSE), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = FALSE, pars.n = NULL, numLoops = 200, numIntervals = 25, optParamListN = 10, optBoundLoops = 10, multicore = FALSE, multicorePackages = NULL, fileTag = NULL, numSimsToAverage = 10,  ...) {

  #split data into different folds
  sampleVec <- rep(1:kFold, ceiling( nrow( data)/kFold ) )
  data$fold <- sample( sampleVec , size=nrow( data) ,  replace=FALSE )

  df.fit.runStats <- NULL
  df.test.runStats <- NULL
  for(fold in 1:kFold) {
    #fit data
    df.fit <- data[data$fold != fold, ]
    #test data
    df.test <- data[data$fold == fold, ]

    #add fold number to the fileTag
    foldNum <- paste("fold-",fold, sep="")
    fileTag.fit <- paste(fileTag,foldNum, "fit", sep="_")
    fileTag.test <- paste(fileTag,foldNum, "test", sep="_")

    #output information into subDir
    curDir <- getwd()
    chutils::ch.newDir (curDir, foldNum)
    foldDir <- getwd()

    ##### fit the model on the fit data
    chutils::ch.newDir (foldDir, "fit")
      run.list <- rrwRunSmartGridSearch(data = df.fit, rrwModelList=rrwModelList, minN = minN, dataOverlapCol = dataOverlapCol, RwSamplesCol = RwSamplesCol, dataRtCol = dataRtCol, correctCol = correctCol, correctVals = correctVals, loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, numLoops = numLoops, numIntervals = numIntervals, optParamListN = optParamListN, optBoundLoops = optBoundLoops, multicore = multicore, multicorePackages = multicorePackages, numSimsToAverage = numSimsToAverage, fileTag = fileTag.fit, ...)

      df.fitted.fit <-run.list$df.fitted
      runStats.fit <- run.list$runStats

      df.fit.stats <- rrwRunStatsToDataframe(runStats.fit)
      df.fit.stats$kFold <- fold
      df.fit.runStats <- chutils::ch.rbind(df.fit.runStats, df.fit.stats)

      rrwPlotSGSoutput(df.fitted.fit, rrwModelList, dataRtCol = "rt", dataPhitCol = "pHit", rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", fileTag = fileTag.fit,numSimsToPlot = numSimsToAverage)
    setwd(foldDir)
    ###### finish fit ######

    ##### test the fit on the test data
    chutils::ch.newDir (foldDir, "test")
      #add fixed parameters to the rrwModelList
      fixedRRWmodelList <- rrwConvertFreeRRWModelToFixedRRWModel(rrwModelList=rrwModelList, runStatsList = runStats.fit)
      #fit the model
      testRun.list <- rrwRunWithFixedParameters(data = df.test, rrwModelList = fixedRRWmodelList, minN = minN, dataOverlapCol = dataOverlapCol, RwSamplesCol = RwSamplesCol, dataRtCol = dataRtCol, correctCol = correctCol, correctVals = correctVals, loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, numSimsToAverage = numSimsToAverage, fileTag = fileTag.test)

      df.fitted.test <-testRun.list$df.fitted
      runStats.test <- testRun.list$runStats
      df.test.stats <- rrwRunStatsToDataframe(runStats.test)
      df.test.stats$kFold <- fold
      df.test.runStats <- chutils::ch.rbind(df.test.runStats, df.test.stats)

      rrwPlotSGSoutput(df.fitted.test, fixedRRWmodelList, dataRtCol = "rt", dataPhitCol = "pHit", rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", fileTag = fileTag.test,numSimsToPlot = numSimsToAverage)
    setwd(foldDir)
    ###### finish test ######

    setwd(curDir)

    #save progress, which gets overwritten as the analysis progresses. It is a safety measure incase of a crash
    filename <- paste(fileTag,"kFoldRunStatsFit.txt")
    write.table(df.fit.runStats, filename, col.names=T, row.names=F, quote=F, sep="\t")
    filename <- paste(fileTag,"kFoldRunStatsTest.txt")
    write.table(df.test.runStats, filename, col.names=T, row.names=F, quote=F, sep="\t")
  }

    meanStat <- data.frame(t(sapply(df.fit.runStats, mean)))
    medianStat <- data.frame(t(sapply(df.fit.runStats, median)))
    sdStat <- data.frame(t(sapply(df.fit.runStats, sd)))

    df.fit.out <- rbind(meanStat, medianStat, sdStat)
    df.fit.out$stat <- c("mean", "median", "sd")
    df.fit.out$kFold <- kFold

    filename <- paste(fileTag,"kFoldRunStatsFit.txt")
    write.table(df.fit.runStats, filename, col.names=T, row.names=F, quote=F, sep="\t")
    filename <- paste(fileTag,"kFoldRunSummaryStatsFit.txt")
    write.table(df.fit.out, filename, col.names=T, row.names=F, quote=F, sep="\t")

    meanStat <- data.frame(t(sapply(df.test.runStats, mean)))
    medianStat <- data.frame(t(sapply(df.test.runStats, median)))
    sdStat <- data.frame(t(sapply(df.test.runStats, sd)))

    df.test.out <- rbind(meanStat, medianStat, sdStat)
    df.test.out$stat <- c("mean", "median", "sd")
    df.test.out$kFold <- kFold

    filename <- paste(fileTag,"kFoldRunStatsTest.txt")
    write.table(df.test.runStats, filename, col.names=T, row.names=F, quote=F, sep="\t")
    filename <- paste(fileTag,"kFoldRunSummaryStatsTest.txt")
    write.table(df.test.out, filename, col.names=T, row.names=F, quote=F, sep="\t")

  out.list <- list(runStatsFit = df.fit.runStats, runStatsTest = df.test.runStats, summaryStatsFit = df.fit.out, summaryStatsTest = df.test.out)
  return(out.list)
}
