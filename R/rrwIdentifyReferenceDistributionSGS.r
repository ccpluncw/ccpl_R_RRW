#' This function takes raw choice data, raw values data, and the rrwModelList as input, and identifies the reference distributions using smartGridSearch, and outputs the smartGridSearch results.
#'
#' Function that takes raw choice data, raw values data, and the rrwModelList as input, and identifies the reference distributions using smartGridSearch, and outputs the smartGridSearch results.
#' @param df.choiceDat A dataframe with the choice data from a typical PVT experiment.
#' @param df.valueDat A dataframe with the raw values data from a typical values experiment.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param refDistGenFn A function to generate the observation for the reference distribution.  This function must take the following arguments: n, mean, sd.  DEFAULT = rnorm.
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
#' @param combFun If there are multiple items contributing to a single distribution, this function describes how the values will be combined across items in both the X and Y distributions. The function must combine elements of a list that might be of different different lengths. The default just flattens the list into one large vector.  DEFAULT = ch.maxAveComb  (with probMax = 0.5)
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = FALSE.
#' @param onlyMinimizeOnPhvoFit A boolean that specifies whether to minimize the function only on the fit of the p(HVO) data (if set to TRUE).  If set to false, the function will be minimized on the basis of both the fit of the p(HVO) and RT data.  DEFAULT = FALSE.
#' @param splitByRefDist A boolean that specifies whether to group the dataset by those trial that are above the reference distribution and those that are below the reference distribution.  DEFAULT = TRUE.
#' @param fitByRefDist A boolean that specifies whether fit the p(HVO) and RT functions by the grouping created when splitByRefDist = TRUE. If splitByRefDist = FALSE, then this argument is irrelevant. Generally, you want splitByRefDist = TRUE and  fitByRefDist = FALSE to find a reference distribution that fits an unbiased observer for a particular condition in your experiment and then see how your manipulations influence either another reference distribution or the response biases in the RRW.  DEFAULT = FALSE.
#' @param splitByCorrect A boolean that specifies whether to use the "correct/incorrect" column as a grouping variable for the RT data when fitting the RT function. When set to FALSE, the RT function will be fit to the collapsed correct and incorrect responses. DEFAULT = TRUE.
#' @param numLoops An integer that specifies the number of loops with different parameter values to be searched in each run of the smartGridSearch before the best fits for that run are identified. DEFAULT = 200.
#' @param numIntervals An integer that specifies the number of intervals that the range between the upper and lower bounds of each parameter space will be segmented into.   DEFAULT = 25.
#' @param optParamListN An integer that specifies the number of best fit runs that will be used to extract the new bounds for next set of grid search loops.  So, if optParamListN == 10, then the parameter values from the 10 best runs will be used to refine the upper and lower bounds of for each parameter. DEFAULT = 10.
#' @param optBoundLoops An integer that specifies the number of times to re-simulate the best fit run in order to get a mean and SD for the minimized fit statistic. If the mean fit statistic is no longer the best fit, the program will move to the next best fit run, and so on until the it identifies the best mean fit run. DEFAULT = 10.
#' @param multicore An boolean that specifies whether the program should be run using multiple processing cores. DEFAULT = FALSE.
#' @param multicorePackages An vector of strings specifying the package names that are used in "fn." This needs to be used when multicore = TRUE. DEFAULT = NULL.
#' @param fileTag A string that is appended to the name of files to identify the analysis and experiment. The default is NULL, whereby the filetag will just be based on a timestamp.
#''
#' @return The file containing the overlaps with the identified reference distribution is saved (<fileTag>Overlaps.txt) and the final parameters of the reference distribution with the fits (<fileTag>FitResults.txt). The function returns a list containing the fitted data to the p(HVO) and RT relative to the overlaps of the reference distribution.
#' @keywords RRW smartGridSearch reference distributions
#' @export
#' @examples rrwIdentifyReferenceDistributionSGS(df.choiceDat=moralsData, df.valueDat = allValues, myModelList, item1cols = "Item1", item2cols = "Item2", respChoiceCol = "keyDef", respChoiceVal = c("Yes", "No"), RTcol = "res.rt",chanceThreshold = 0.5,  lowRTquantileThreshold = 0.025, highRTquantileThreshold = 0.975, minOverlapN = 20, pars.n = 3, equalizeRTandPhit = TRUE, minimizeStat = "BIC", roundThreshold = 0.1, roundDirection = ceiling, overlapNumRuns = 1000, combFun = ch.maxAveComb, numLoops = 2500, numIntervals = 100, optParamListN = 10, optBoundLoops = 10, multicore = TRUE, multicorePackages = c('RRW'), fileTag = NULL, probMax = 0.5)

rrwIdentifyReferenceDistributionSGS <- function(choiceData, valueData, rrwModelList, refDistGenFn = rnorm, item1cols, item2cols, respChoiceCol, respChoiceVal, RTcol, chanceThreshold = 0, lowRTquantileThreshold = 0, highRTquantileThreshold = 0, minOverlapN = 0, pars.n = NULL, equalizeRTandPhit = TRUE, minimizeStat = "BIC", roundThreshold = 0.1, roundDirection = ceiling, numOverlapBins = 10, overlapNumRuns = 1000, combFun = ch.maxAveComb, useTwoParameterModel = FALSE, onlyMinimizeOnPhvoFit = FALSE, splitByRefDist = TRUE, fitByRefDist = FALSE, splitByCorrect = TRUE, numLoops = 200, numIntervals = 25, optParamListN = 10, optBoundLoops = 10, multicore = FALSE, multicorePackages = NULL, fileTag = NULL, ...) {

  #if fileTag does not exist, create one based on the timestamp
  if(is.null(fileTag)) {
    fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), "RRW", sep="_")
  }

  #get parameter information
  pList <- rrwGetModelParameterSpecs(rrwModelList)

  #add the necessary parameters to the "otherParameterList"
  otherParamList <- list(refDistGenFn = refDistGenFn, item1cols = item1cols, item2cols = item2cols, respChoiceCol = respChoiceCol, respChoiceVal = respChoiceVal, RTcol = RTcol, chanceThreshold = chanceThreshold, lowRTquantileThreshold=lowRTquantileThreshold, highRTquantileThreshold=highRTquantileThreshold, minOverlapN=minOverlapN, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, roundThreshold=roundThreshold, roundDirection=roundDirection, numOverlapBins=numOverlapBins, overlapNumRuns=overlapNumRuns, overlapBootstrapMulticore = FALSE, useTwoParameterModel = useTwoParameterModel, onlyMinimizeOnPhvoFit = onlyMinimizeOnPhvoFit, splitByRefDist = splitByRefDist, fitByRefDist = fitByRefDist, splitByCorrect = splitByCorrect, combFun=combFun, ...)

  #move the bounds, columns, and pars.n from pList to the appropriate format
  pU <- pList$pU
  pL <- pList$pL
  pI <- pList$pI
  if(is.null(pars.n)) {
    otherParamList <- c(otherParamList, pList$pars.n)
  } else {
    otherParamList <- c(otherParamList, pars.n = pars.n)
  }
  otherParamList <- c(otherParamList, list(df.choiceDat = choiceData, df.valueDat = valueData))

  #get the start time
  start <- Sys.time()

  #run the smartGridSearch to identify the optimal parameter values
  #xgrid gives the best fit parameters
  x.grid <- smartGridSearch::smartGridSearch(rrwAssessComparisonDistribitionFit, pU, pL, pI, otherParamList, numLoops = numLoops, numIntervals = numIntervals, optParamListN = optParamListN, optBoundLoops = optBoundLoops, multicore = multicore, multicorePackages = multicorePackages)

  ## get the end time
  end <- Sys.time()
  # see how long the analysis took to run
  diff <- end-start
  #output the effects
  filename <- paste(fileTag,"GridOut.txt")
  sink(filename)
    print(x.grid)
    print(diff)
  sink(NULL)


  #create an input list to run "getMeanRRWfit" tht will run the simulation with the best fit
  #parameters multiple times and calculate the average of the runs.  The average is important to
  #smooth out the random fluctuations of the stochastic process.
  x.in <- as.list(x.grid$final)
  #NULL out the unused parameters
  x.in$r_out <- NULL
  x.in$sdR_out <- NULL
  otherParamList$minimizeStat <- NULL
  #add the output filenames
  x.in$sinkFilename <- paste(fileTag,"FitResults.txt")
  #add the other paramters
  x.in <- c(x.in, otherParamList)
  #How many simulations do you run and average together?
  x.in$overlapOutFile <- paste(fileTag,"Overlaps.txt")
  x.in$figureOutFile <- paste(fileTag,"RTpHVOfits.pdf")
  x.in$overlapBootstrapMulticore <- TRUE
  x.in$overlapNumRuns <- x.in$overlapNumRuns*10

  #save the Parameter values as an r data format
  filename <- paste(fileTag,"x.in.RData")
  save(x.in, file=filename)
  #run the RRW to get the fit given the parameters output by the smartGridSearch
  df.fitted <- do.call(getReferenceDistributionFit, x.in)

  #add the minimization statistic to the output file
  sink(x.in$sinkFilename, append = TRUE)
    cat("\n Minimization Statistic = ", minimizeStat, "\n")
  sink(NULL)

 return(df.fitted)
}
