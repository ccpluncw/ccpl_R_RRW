#' A function to assess how well a specific comparison distribution fit with choice data
#'
#' This function assesses how well a specific comparison distribution fit with choice data. The comparison distribution is normal, with a mean and SD.
#' @param df.choiceDat A dataframe with the choice data from a typical PVT experiment.
#' @param df.valueDat A dataframe with the raw values data from a typical values experiment.
#' @param rdN An integer specifying the number of observations in the reference distribution.
#' @param rdMean An number specifying the mean of the reference distribution (Gaussian).
#' @param refDistGenFn A function to generate the observation for the reference distribution.  This function must take the following arguments: n, mean, sd.  DEFAULT = rnorm.
#' @param rdSD An number specifying the standard deviation of the reference distribution (Gaussian).
#' @param item1cols a vector of strings that specifies the names of the columns in "df.choiceDat" that contains the the probes in Item 1.
#' @param item2cols a vector of strings that specifies the names of the columns in "df.choiceDat" that contains the the probes in Item 2.
#' @param respChoiceCol a string that specifies the name of the column in "df.choiceDat" that contains the the participant's response to the prompt - yes take action or no take no action.
#' @param respChoiceVal a vector of two values that specifies the choose Item1 option ("yes" take action in many morals experimants, thus saving Item1) value (index 1) and the choose Item1 option ("no" take no action in many morals experimants, thus saving Item2) value (index 2).
#' @param RTcol a string that specifies the name of the column in "df.choiceDat" that contains the RT for each trial.
#' @param chanceThreshold A number specifying the quantile that all subjects whose p(HVO) falls below it will be removed from the dataset. DEFAULT = 0 (none removed)
#' @param lowRTquantileThreshold A number specifying the quantile that, for each overlapRound, individual RTs that fall below it will be removed from the dataset. DEFAULT = 0.0 (none removed)
#' @param highRTquantileThreshold A number specifying the quantile that, for each overlapRound, individual RTs that fall above  it will be removed from the dataset. DEFAULT = 1.0 (none removed)
#' @param minOverlapN A number specifying the minimum number of responses an overlapRound condition must have to remain in the dataset. DEFAULT = 0 (none removed)
#' @param pars.n The number of free parameters. When NULL, the program will attempt to calculate the number of free parameters from the input. Default = 3.
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.
#' @param minimizeStat A string that specifies which statistic to minimize when optimizing the model fit.  The options are: "BIC" , "AIC" , or "R_Square". Default is "BIC".
#' @param roundThreshold An integer that specifies the nearest interval that the overlaps should be rounded to. DEFAULT = 0.1 (round to the nearest 0.1)
#' @param roundDirection An option that specifies the rounding direction: ceiling (always round up), floor (always round down), or round (round to the nearest value, up or down). DEFAULT = ceiling
#' @param numOverlapBins An integer that specifies the number of bins the overlaps should be binned into. This is only used if roundThreshold = NULL.  DEFAULT = 10 (round overlaps so that there are 10 bins)
#' @param overlapNumRuns the number of runs to do in the bootstrap that calculates the overlap of the reference distribution and the choice itms. DEFAULT = 1000.
#' @param overlapBootstrapMulticore A boolean that specifies whether to run the bootstrap that calculates the value overlaps in multicore mode (in parallel).  This should be set to FALSE if you are using this function within the smartGridSearch. DEFAULT = FALSE.
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = FALSE.
#' @param onlyMinimizeOnPhvoFit A boolean that specifies whether to minimize the function only on the fit of the p(HVO) data (if set to TRUE).  If set to false, the function will be minimized on the basis of both the fit of the p(HVO) and RT data.  DEFAULT = FALSE.
#' @param splitByRefDist A boolean that specifies whether to group the dataset by those trial that are above the reference distribution and those that are below the reference distribution.  DEFAULT = TRUE.
#' @param fitByRefDist A boolean that specifies whether fit the p(HVO) and RT functions by the grouping created when splitByRefDist = TRUE. If splitByRefDist = FALSE, then this argument is irrelevant. Generally, you want splitByRefDist = TRUE and  fitByRefDist = FALSE to find a reference distribution that fits an unbiased observer for a particular condition in your experiment and then see how your manipulations influence either another reference distribution or the response biases in the RRW.  DEFAULT = FALSE.
#' @param splitByCorrect A boolean that specifies whether to use the "correct/incorrect" column as a grouping variable for the RT data when fitting the RT function. When set to FALSE, the RT function will be fit to the collapsed correct and incorrect responses. DEFAULT = TRUE.
#' @param combFun If there are multiple items contributing to a single distribution, this function describes how the values will be combined across items in both the X and Y distributions. The function must combine elements of a list that might be of different different lengths. The default just flattens the list into one large vector.  DEFAULT = ch.maxAveComb  (with probMax = 0.5)
#' @keywords assess reference distribution fit
#' @return The minimization statistic for the fit of the model to the RT and pHit (proportion crossed each boundary) data.  This is the value that will be miniized by the optimization program.
#' @export
#' @examples rrwAssessComparisonDistribitionFit (df.choiceDat=moralsData, df.valueDat = allValues, rdN=100, rdMean=3, rdSD=1, item1cols = "Item1", item2cols = "Item2", respChoiceCol = "keyDef", respChoiceVal = c("Yes", "No"), RTcol = "res.rt",chanceThreshold = 0.5,  lowRTquantileThreshold = 0.025, highRTquantileThreshold = 0.975, minOverlapN = 20, pars.n = 3, equalizeRTandPhit = TRUE, minimizeStat = "BIC", roundThreshold = 0.1, roundDirection = ceiling, numOverlapBins = 10, overlapNumRuns = 1000, combFun = ch.maxAveComb, probMax = 0.5)

rrwAssessComparisonDistribitionFit <- function (df.choiceDat, df.valueDat, rdN, rdMean, rdSD, refDistGenFn = rnorm, item1cols, item2cols, respChoiceCol, respChoiceVal, RTcol, chanceThreshold = 0, lowRTquantileThreshold = 0, highRTquantileThreshold = 0, minOverlapN = 0, pars.n = 3, equalizeRTandPhit = TRUE, minimizeStat = "BIC", roundThreshold = 0.1, roundDirection = ceiling, numOverlapBins = 10, overlapNumRuns = 1000, overlapBootstrapMulticore = FALSE, useTwoParameterModel = FALSE, onlyMinimizeOnPhvoFit = FALSE, splitByRefDist = TRUE, fitByRefDist = FALSE, splitByCorrect = TRUE, combFun = ch.maxAveComb, ...) {

  #make sure minimizeStat is valid
  minimizeOpts <- c("BIC", "AIC", "R_Square")
  if(!(minimizeStat %in% minimizeOpts) ) {
    stop (paste("you set minimizeStat to:", minimizeStat, ", but it must be one of the following:", minimizeOpts, sep=" "))
  }

  #check if parameters are valid
  if( (rdN > 0) & (rdSD >= 0)) {
    #generate reference distribution.
    refDist <- refDistGenFn (n = rdN, mean = rdMean, sd = rdSD)

    df.refDist <- data.frame(prompt = "refDist", respS = refDist)

    #paste them together
    df.valueAll <- rbind(df.valueDat, df.refDist)

    #create a column in df.prompts for the reference distribution
    df.prompts <- chMorals::ch.getPrompts(df.choiceDat, item1cols = item1cols, item2cols = item2cols)
    #calculate the overlaps relative to the reference distribution
    df.over <- chValues::ch.MCbatchOverlapPromptFile(df.valueAll$respS, df.valueAll$prompt, df.prompts = df.prompts, itemAcolNames = item1cols, itemBcolNames = item2cols, numRuns = overlapNumRuns, outFile = NULL,  multicore = overlapBootstrapMulticore, combFun = combFun, ...)

    overlapItem1cols <- paste("IA", seq(1,length(item1cols),1), sep="")
    overlapItem2cols <- paste("IB", seq(1,length(item2cols),1), sep="")

    #add overlaps
    df.out1 <- chMorals::ch.mergeChoiceDataWithOverlapsData(df.choiceDat, df.over, "overlap", "direction", respChoiceCol = respChoiceCol, respChoiceVal = respChoiceVal, item1cols = item1cols, item2cols = item2cols, overlapItem1cols = overlapItem1cols, overlapItem2cols = overlapItem2cols, roundThreshold = roundThreshold, roundDirection = roundDirection, numOverlapBins= numOverlapBins,outfile = NULL)

    #filter with overlaps
    df.out2 <- chMorals::ch.filterOverlapsByQuantile(df.out1, "sn", RTcol = RTcol, overlapRoundCol = "overlapRound", correctCol = "correct", correctVals = c(1,0), chanceThreshold = chanceThreshold, lowRTquantileThreshold = lowRTquantileThreshold, highRTquantileThreshold = highRTquantileThreshold, minOverlapN = minOverlapN, statsOutputFile = NULL)

    #create grouping variable whereby the reference distribution is either the HVO or LVO
    grpColSplit <- NULL
  	grpColsFit <- NULL
  	if(splitByRefDist) {
      df.out2$refValue<-ifelse(df.out2$dirOverlap <= 0, "refHVO", "refLVO")
  		grpColSplit <- "refValue"
      if(fitByRefDist) {
    		grpColsFit <- "refValue"
    	}
  	}

  	df.out3 <- rrwRawDataToRRWFormatData(df.out2, grpColSplit, dataRtCol = RTcol, correctCol = "correct", correctVals = c(1,0), dataOverlapCol= "overlapRound", splitByCorrect = splitByCorrect)

    if(splitByCorrect) {
      correctCol = "correct"
      correctVals = c(TRUE, FALSE)
    } else {
      correctCol = NULL
      correctVals = NULL
    }

    regFits <- getPredictedRTpHVOfitForRefDist(df.out3, grpCols = NULL,RTcol = "rt", pHVOcol = "pHit", overlapRoundCol="overlapRound", refDistCol = grpColsFit, correctCol = correctCol, correctVals = correctVals, nCol = "n", minNperOverlap = minOverlapN, useTwoParameterModel = useTwoParameterModel)

    #get the fit statistics
    if(onlyMinimizeOnPhvoFit) {
      out.rss <- getMinR2(as.vector(regFits$data$pHVOfit), as.vector(regFits$data$pHit), standardize = FALSE)
      out.BIC <- chutils::ch.BIC(as.vector(regFits$data$pHit), as.vector(regFits$data$pHVOfit), pars.n, standardize = FALSE)
      out.AIC <- chutils::ch.AIC(as.vector(regFits$data$pHit), as.vector(regFits$data$pHVOfit), pars.n, standardize = FALSE)

    } else {
      if(!is.null(correctCol)) {
        pHVOdata <- regFits$data[regFits$data[[correctCol]] == correctVals[1], ]
      } else {
        pHVOdata <- regFits$data
      }
      out.rss <- getMinR2RRW(as.vector(pHVOdata$pHVOfit), as.vector(pHVOdata$pHit), as.vector(regFits$data$RTfit), as.vector(regFits$data$rt), equalizeRTandPhit = equalizeRTandPhit)
      out.BIC <- getICRRW(as.vector(pHVOdata$pHVOfit), as.vector(pHVOdata$pHit), as.vector(regFits$data$RTfit), as.vector(regFits$data$rt), pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "BIC")
      out.AIC <- getICRRW(as.vector(pHVOdata$pHVOfit), as.vector(pHVOdata$pHit), as.vector(regFits$data$RTfit), as.vector(regFits$data$rt), pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "AIC")
    }

#This sink is for debugging
# sink("debug1.txt", append = T)
#   cat("\n\nrdM: ", rdMean, "; mean(refDist): ", mean(refDist), " ; rdSD: ", rdSD, " ; sd(refDist): ", sd(refDist), "\n")
#   cat("\n\n ******** p(HVO) fit ******** \n\n")
#   for(i in 1:length(regFits$pHVOFit)) {
#     print(summary(regFits$pHVOFit[[i]]$nlsObject))
#   }
#   cat("\n\n ******** p(HVO) R Square ******** \n\n")
#   for(i in 1:length(regFits$pHVOFit)) {
#     print(regFits$pHVOFit[[i]]$r2)
#   }
#   cat("\n\n ******** RT fits ******** \n\n")
#   for(i in 1:length(regFits$RTfit)) {
#     print(summary(regFits$RTfit[[i]]))
#   }
# sink(NULL)

  } else {
    #if the parameters are not valid, then assign the minimization variable Infinite value (Inf)
    out.rss <- Inf
    out.BIC <- Inf
    out.AIC <- Inf
  }

  switch(minimizeStat, BIC = return(out.BIC), AIC = return(out.AIC), R_Square = return(out.rss))


}
