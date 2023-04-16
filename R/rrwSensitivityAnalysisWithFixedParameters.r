#' This function runs a sensisity analysis with a set of fixed parameters (specified in a rrwModelList) to generate the fits to the data. It takes the raw data and the rrwModelList as input, collapses the data, codes the columns, gets the parameter information, runs the getMeanRRWfit, and outputs the smartGridSearch results. rrwModelList should only have fixed parameter values (upperBound==lowerBound). The sensitivity analysis will randomly drop probes from the probe list and check the fit.
#'
#' Function that runs a set of fixed parameters (specified in a rrwModelList) to generate the fits to the data.
#' @param data This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param probeCols a vector of strings that specifies the names of the columns in "data" that contains the the probes.
#' @param numProbesToRemove A vector of integers that specifies the number of probes to remove with each run. DEFAULT = c(1); which is a "leave one out" protocol.
#' @param maxNumCombinations An integer that specifies the maximum number of combinations to remove each run. When the eniter set is less than maxNumCombinations, the entire set of combinantions is removed. But often when the value of numProbesToRemove is large, the number of possible combinations gets very large, so it is best to randomly choose a subset of those to run. maxNumCombinations provides an upper limit on the random number of combinations. aDEFAULT = 200
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
#' @param savePlots A boolean specifying whether to save plots for each sensitivity run. Default = FALSE.
#''
#' @return Several files are saved, including <fileTag>rrwData.txt (which contains the inputData), <fileTag>FitResults.txt (which contains the final parameter fits), <fileTag>Fitted.txt (which contains the input data + simulation fit data). The function returns a list containing (1) a list of the parameter values and fit statistics [runStats], and (2) a dataframe containing the input data with the best fit prediction based on the parameter values [df.fitted].
#' @keywords RRW fixed parameter run getMeanRRWfit
#' @export
#' @examples rrwRunWithFixedParameters(tmp.df, tmpModelList, minN = 40, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct01", correctVals = c(1,0), loopsPerRWstep = 200, minimizeStat = 'BIC', equalizeRTandPhit = TRUE, numSimsToAverage = 40, fileTag = NULL)

rrwSensitivityAnalysisWithFixedParameters <- function(data, rrwModelList, probeCols = NULL, numProbesToRemove = c(1), maxNumCombinations = 200, minN = 20, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "rt", correctCol = "correct",correctVals = c(TRUE, FALSE), loopsPerRWstep = 200, minimizeStat = 'R_Square', equalizeRTandPhit = FALSE, pars.n = NULL, fileTag = NULL, numSimsToAverage = 10, savePlots = FALSE, ...) {


    if(is.null(probeCols)) {
      stop("probeCols must not be NULL (it must be a vector of strings)")
    }

    #get the probes
    probes <- vector()
    for(i in probeCols) {
      probes <- append(probes, unique(data[,i]))
    }
    probes <- unique(probes)
    numProbes <- length(probes)
    numProbeCols <- length(probeCols)

    df.runStats <- NULL
    for(rn in numProbesToRemove) {
      cat("run ", rn, "of ", numProbesToRemove, " probes to remove \n")
      #get all the combinations of probes for this value of numProbesToRemove
      runItemVec <- chValues::ch.combnVector (numProbes, rn, rn, F)
      #get the total number of combincations for this value of numProbesToRemove
      numItemVecRows <- nrow(runItemVec)

      #randomly sample from the set of items if there are more than maxNumCombinations
      if(numItemVecRows > maxNumCombinations) {
        ind <- seq(1:numItemVecRows)
        ind.samp <- sample(ind, maxNumCombinations, replace = F)
        runItemVec <- runItemVec[ind.samp,]
        numItemVecRows <- nrow(runItemVec)
      }

      sensRN <- paste("sense",rn, sep="-")
      curDir <- getwd()
      chutils::ch.newDir (curDir, sensRN)
      sensDir <- getwd()

      #for each combination
      for(niv in 1:numItemVecRows) {
        cat("iteration ", niv, "of ", numItemVecRows, " interations for ", rn, " probes to remove \n")
        sensNum <- paste(sensRN,niv, sep="-")
        fileTag.tmp <- paste(fileTag,sensNum, sep="_")

        #identify the probes to remove
        removeProbes <- probes[unlist(runItemVec[niv,])]

        ### remove the probes
        df.tmp <- data
        for(pc in probeCols) {
          df.tmp <- df.tmp[!(df.tmp[[pc]] %in% removeProbes ),]
        }

        #run the analysis on the data with the probes removed
        run.list <- rrwRunWithFixedParameters(data = df.tmp, rrwModelList = rrwModelList, minN = minN, dataOverlapCol = dataOverlapCol, RwSamplesCol = RwSamplesCol, dataRtCol = dataRtCol, correctCol = correctCol,correctVals = correctVals, loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, pars.n = pars.n, fileTag = fileTag.tmp, numSimsToAverage = numSimsToAverage, ...)

        df.fitted <-run.list$df.fitted

        if(savePlots) {
          rrwPlotSGSoutput(df.fitted, rrwModelList, dataRtCol = "rt", dataPhitCol = "pHit", rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", fileTag = fileTag.tmp, numSimsToPlot = numSimsToAverage)
        }

        runStats <-run.list$runStats

        df.tmp.stats <- rrwRunStatsToDataframe(runStats)
        df.tmp.stats$numProbesRemoved <- rn
        df.tmp.stats$probesRemoved <- paste(removeProbes,collapse=";")
        df.runStats <- chutils::ch.rbind(df.runStats, df.tmp.stats)

      }
      setwd(curDir)

      #save progress, which gets overwritten as the analysis progresses. It is a safety measure incase of a crash
      filename <- paste(fileTag,"sensitivityAnalysisRunStats.txt")
      write.table(df.runStats, filename, col.names=T, row.names=F, quote=F, sep="\t")
    }

    df.summary <- data.frame(df.runStats %>% group_by (numProbesRemoved) %>% summarize (mAIC = mean(AIC, na.rm = T),sdAIC = sd(AIC, na.rm = T), mBIC = mean(BIC, na.rm = T), sdBIC = sd(BIC, na.rm = T), mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T)))
    filename <- paste(fileTag,"sensitivityAnalysisRunStats.txt")
    write.table(df.runStats, filename, col.names=T, row.names=F, quote=F, sep="\t")
    filename <- paste(fileTag,"sensitivityAnalysisSummaryStats.txt")
    write.table(df.summary, filename, col.names=T, row.names=F, quote=F, sep="\t")

    out.list <- list(df.runStats = df.runStats, df.summary = df.summary)
    return(out.list)
  }
