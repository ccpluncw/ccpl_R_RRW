#' This function generates an RRW simulation, fits it to empirical data, and and outputs the fit.
#'
#' Function that generates an RRW simulation, fits it to empirical data, and outputs the fit.  It is generally run after the assessRRWfit() is fit with an optimation routine (e.g., smartGridSearch()) to identify the optimal parameters.  The optimation routine should output the parameters that generate the best fit, and those parameters should be input into plotRRWfit().  This takes the same parameters as assessRRWfit().
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param b A vector of number(s) specifying the boundary distance from a 0 startpoint. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bCols" argument. b is specific to the RRW simulation and has no default value
#' @param bcs A vector of positive values that specify the likelihood that the boundary will be reduced because there it is taking too many samples before reaching a threshold. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bcsCols" argument. Essentially, this models a decision-maker reducing the amount of information that they need before a response because it is taking to long.  The boundary reduces when there is little change in the likelihood of a response after N number of steps.  The routine makes this determination as a function of the size of the unaltered boundary (boundaryChangeSensitivity * boundary). Small values (e.g., 0.1), result in an impatient person (relatively early chnage in the boundary values).  Large values (e.g., 0.9), result in a patient person (a boundary that is relatively resistent to change.  A value of 0 is an exception: it will result in no boundary reduction under any circumstances. Impatience increases the likelihood of an error, so it increases the average time of error responses.  Default is 0.25.
#' @param Ter A vector of positive values that adjusts the time of encoding/response (non-decision time). If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter. All effect/dummy codes must >=0. The column names must be specified in the "TerCols" argument. When Ter = 0, a single, best fitting Ter will be calculated. Use this parameter when you have several conditions and you hypothesize that some conditions should have a larger Ter than others.  For example, this may model longer encoding process when different numbers of items are on the screen for different conditions. The hypothetical fastest condition should have Ter = 0, which simply indicates that the RRW program will find the best fit Ter. The slower conditions will have larger Ter values and will be fit ralative to the fastest condition. The change in samples resulting from the Ter effect is a function of the size of the boundary (Ter * boundary). Default = 0, This is still experimental.
#' @param s A vector of signed number(s) between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.   If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of s. The column names must be specified in the "sCols" argument. s = 0 is the default and represents an unbiased start point.
#' @param nSD A vector of positive number(s) representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step.  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD. The column names must be specified in the "nSDCols" argument. nSD = 0 is the default and represents no noise being added to each step.
#' @param db A vector of signed number(s) representing the beta value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of db. The column names must be specified in the "dbCols" argument. db = 0 is the default and represents no IAB.
#' @param da A vector of positive number(s) representing the asymptote value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of da. The column names must be specified in the "daCols" argument. da = 0.2 is the default.
#' @param vc A vector of signed number(s) that indicates the change of overlap that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of vc. The column names must be specified in the "vcCols" argument. The default is 0 because there is no valueChangeEffect.
#' @param bCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of boundary (b).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of b.
#' @param bcsCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the boundaryChangeSensitivity (bcs).  These columns must be effect coded, because the values contained in this column will be multiplied by the value of bcs. All effect/dummy codes must be >=0.
#' @param TerCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the non-decision time (Ter).  These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter.  All effect/dummy codes must be >=0.
#' @param sCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of startValue (s).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of s.
#' @param nSDCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of noiseSD (nSD).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD.
#' @param dbCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decayBeta (db).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of db.
#' @param daCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decay asymptote (da).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of da.
#' @param vcCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of value change (vc).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of overlap.
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param RwSamplesCol A string that identifies the name of the column in the RRW simulations that contains the summary of the samples that you want to use as a simulation for RT.  The possible columns are: "Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile).  The default is "Q50" because the median is more robust than the mean.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param dataCorrectCol A string that identifies the name of the column in data that identifies whether the trials were correct (TRUE) or incorrect (FALSE). The default is "correct"
#' @param loopsPerRWstep A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @param sinkFilename A string that identifies the name of file (.txt) in which the fit statistics will be saved. The default is NULL, whereby the fit statistics will not be saved.
#' @param pars.n The number of free parameters in the model. The default is NULL, whereby the function will calculate the number of free parameters automatically.
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.
#''
#' @return The function returns a list containing (1) a list of the parameter values and fit statistics [runStats], and (2) a dataframe [df.fitted] that contains the "data" plus the fitted values from the model ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the fitted pHit from the model - probability of crossing each boundary); rtFit (the fitted rt values from the model))
#' @keywords RRW random walk plot output fit
#' @export
#' @examples getRRWfit (data, b=14, s=0.1, loopsPerRWstep = 2000, sinkFilename = "outStats.txt")

getRRWfit <- function(data, b, bcs = 0.25, Ter = 0, s=0, nSD=0, db=0, da=0.2, vc = 0, bCols = NULL, bcsCols = NULL, TerCols = NULL, sCols = NULL, nSDCols = NULL, dbCols = NULL, daCols = NULL, vcCols = NULL, dataOverlapCol = "overlap", RwSamplesCol = "Q50", dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = 2000, sinkFilename = NULL, pars.n = NULL, equalizeRTandPhit = FALSE) {

  #the RW columns are defined by the program. The column used to compare with the data RT is a choice of the user.
  RWkeepColumns <- c("overlap", RwSamplesCol, "pCross", "correct")
  mergeByDataColumns <- c(dataOverlapCol, dataCorrectCol)
  mergeByRWColumns <- c("overlap", "correct")

  #add the parameter effect columns. if they are NULL, they don't appear anyways
  RWkeepColumns <- c(RWkeepColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols)
  mergeByDataColumns <- c(mergeByDataColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols)
  mergeByRWColumns <- c(mergeByRWColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols)

  df.fitted <- getPredictedRRWpoints(data = data, RWkeepColumns = RWkeepColumns, mergeByDataColumns = mergeByDataColumns, mergeByRWColumns = mergeByRWColumns, dataRtCol = dataRtCol, RwSamplesCol = RwSamplesCol, dataOverlapCol = dataOverlapCol, b = b, boundaryChangeSensitivity = bcs, Ter = Ter, startValue = s,  noiseSD = nSD, decayBeta = db, decayAsymptote = da, valueChange = vc, bCols = bCols, bcsCols = bcsCols, TerCols = TerCols, sCols = sCols, nSDCols = nSDCols, dbCols = dbCols, daCols = daCols, vcCols = vcCols, loops = loopsPerRWstep)

  Q50.lm <- lm(df.fitted[[dataRtCol]]~df.fitted[[RwSamplesCol]])

  #caluculate (1-r2)
  #get the  minimization variable (1-r2) for the pHit
  pCor.rss <- getMinR2(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol])

  #and RT.
  Q50.rss <- getMinR2(df.fitted$rtFit,df.fitted[[dataRtCol]])

  #combine (1-r2) for pHit and RT
  ave.rss <- (Q50.rss + pCor.rss)/2

  out.rss <- getMinR2RRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], equalizeRTandPhit = equalizeRTandPhit)

  #calculate BIC
  if(is.null(pars.n)) {
    pars.n = 1 + length(b) + ifelse(length(bcs)==1 & bcs[1]==0.25, 0, length(bcs)) + ifelse(length(Ter)==1 & Ter[1]==0, 0, length(Ter)) + ifelse(length(s)==1 & s[1]==0, 0, length(s)) + ifelse(length(nSD)==1 & nSD[1]==0, 0, length(nSD)) + ifelse(length(db)[1]==1 & db==0, 0, length(db)) + ifelse(length(da)[1]==1 & da==0.2, 0, length(da)) + ifelse(length(vc)==1 & vc[1]==0, 0, length(vc))
  }
  out.BIC <- getICRRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "BIC")
  out.AIC <- getICRRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "AIC")

  if (!is.null(sinkFilename)) {
    sink(sinkFilename)
      cat("\n\n **************** RRW Statistics **************** \n\n")
      cat("\n\n ******** Parameter Values ******** \n\n")

      cat(" Non-decision time (Ter) for the condition when the dummy/effect codes of Ter=0 (scaled to RT).\n")
      cat(" RTi (Ter 0) = ", coef(Q50.lm)[1], "\n\n")

      cat(" Non-decision time (Ter) for the condition(s) when the dummy/effect codes of Ter != 0 (scaled as: Ter * boundary):\n")
      if(!is.null(TerCols)) cat(" Ter Effect Codes = ", TerCols, "\n")
      cat(" Ter = ", Ter, "\n\n")

      if(!is.null(bCols)) cat(" Boundary Effect Codes = ", bCols, "\n")
      cat(" Boundary = ", b, "\n\n")

      if(!is.null(bcsCols)) cat(" BoundaryChangeSensitivity Effect Codes = ", bcsCols, "\n")
      cat(" BoundaryChangeSensitivity = ", bcs, "\n\n")

      if(!is.null(sCols)) cat(" StartValue Effect Codes = ", sCols, "\n")
      cat(" StartValue = ", s, "\n\n")

      if(!is.null(nSDCols)) cat(" NoiseSD Effect Codes = ", nSDCols, "\n")
      cat(" NoiseSD = ", nSD, "\n\n")

      if(!is.null(dbCols)) cat(" DecayBeta Effect Codes = ", dbCols, "\n")
      cat(" DecayBeta = ", db, "\n\n")

      if(!is.null(daCols)) cat(" DecayAsymptote Effect Codes = ", daCols, "\n")
      cat(" DecayAsymptote = ", da, "\n\n")

      if(!is.null(vcCols)) cat(" ValueChange Effect Codes = ", vcCols, "\n")
      cat(" ValueChange = ", vc, "\n\n")


      cat("\n\n ******** Model Fit Statistics ******** \n\n")
      cat(" The Model fit statistics are based on the simultaneous fit of both the RT and p(HVO). \n")
      cat(" To assess the final fit, Equalize pHit and RT = ", equalizeRTandPhit, "\n\n")

      cat(" N Parameters = ", pars.n, "\n\n")

      cat(" Model R Square = ", 1 - out.rss, "\n")
      cat(" Model BIC = ", out.BIC, "\n")
      cat(" Model AIC = ", out.AIC, "\n")

      cat("\n\n **************** For informational use only **************** \n\n")
      cat(" Below are the R squares for the Model's fit to p(HVO) and RT separately.\n\n")

      cat(" p(HVO) R Square = ", 1 - pCor.rss, "\n")
      cat(" RT R Square = ", 1 - Q50.rss, "\n\n")

      cat(" Below is the Simulation Epoch Transformation. It scales the number of random walk samples to RT.\n")
      cat(" RTb = ", coef(Q50.lm)[2], "\n\n")

    sink(NULL)
  }

  fitStats <- list(AIC = out.AIC, BIC = out.BIC, r2 = 1 - out.rss)
  ter_0.val <- coef(Q50.lm)[1]
  if(!is.null(bCols)) {
    boundary <- list(columns = bCols, values = b)
  } else {
    boundary <- list(columns = "b", values = b)
  }
  if(!is.null(bcsCols)) {
    boundaryChangeSensitivity <- list(columns = bcsCols, values = bcs)
  } else {
    boundaryChangeSensitivity <- list(columns = "bcs", values = bcs)
  }
  if(!is.null(TerCols)) {
    Ter <- list(columns = c("ter_0", TerCols), values = c(ter_0.val, Ter))
  } else {
    Ter <- list(columns = "ter_0", values = ter_0.val)
  }

  if(!is.null(sCols)) {
    StartValue <- list(columns = sCols, values = s)
  } else {
    StartValue <- list(columns = "s", values = s)
  }
  if(!is.null(nSDCols)) {
    NoiseSD <- list(columns = nSDCols, values = nSD)
  } else {
    NoiseSD <- list(columns = "nSD", values = nSD)
  }
  if(!is.null(dbCols)) {
    DecayBeta <- list(columns = dbCols, values = db)
  } else {
    DecayBeta <- list(columns = "db", values = db)
  }
  if(!is.null(daCols)) {
    DecayAsymptote <- list(columns = daCols, values = da)
  } else {
    DecayAsymptote <- list(columns = "da", values = da)
  }
  if(!is.null(vcCols)) {
    ValueChange <- list(columns = vcCols, values = vc)
  } else {
    ValueChange <- list(columns = "vc", values = vc)
  }
  runInfo <- list(freeParameters = pars.n, numSimRuns = 1, equalizePhitRT = equalizeRTandPhit)
  parameters <- list(Ter = Ter, b = boundary, bcs = boundaryChangeSensitivity, s = StartValue, nSD = NoiseSD, db = DecayBeta, da = DecayAsymptote, vc = ValueChange)
  runStats <- list(parameters = parameters, fitStats = fitStats, runInfo = runInfo)

  outlist <- list(df.fitted = df.fitted, runStats = runStats)


  return(outlist)
}
