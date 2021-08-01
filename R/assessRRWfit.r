#' This function fits the RRW and returns (1-r2).
#'
#' Function that fits the RRW and returns (1-r2), where r2 is the fit of the RRW simulation to the empirical RT and error data
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param b A vector of number(s) specifying the boundary distance from a 0 startpoint. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bCols" argument. b is specific to the RRW simulation and has no default value
#' @param s A vector of signed number(s) between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.   If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of s. The column names must be specified in the "sCols" argument. s = 0 is the default and represents an unbiased start point.
#' @param nSD A vector of positive number(s) representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step.  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD. The column names must be specified in the "nSDCols" argument. nSD = 0 is the default and represents no noise being added to each step.
#' @param db A vector of signed number(s) representing the beta value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of db. The column names must be specified in the "dbCols" argument. db = 0 is the default and represents no IAB.
#' @param da A vector of positive number(s) representing the asymptote value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of da. The column names must be specified in the "daCols" argument. da = 0.2 is the default.
#' @param vc A vector of signed number(s) that indicates the change of overlap that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of vc. The column names must be specified in the "vcCols" argument. The default is 0 because there is no valueChangeEffect.
#' @param bCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of boundary (b).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of b.
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
#' @param minimizeStat A string that specifies which statistic to minimize when optimizing the model fit.  The options are: "BIC" , "AIC" , or "R_Square". Default is "BIC".
#' @param pars.n The number of free parameters. When NULL, the program will attempt to calculate the number of free parameters from the input. Default = NULL.
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.
#''
#' @return The minimization statistic for the fit of the model to the RT and pHit (proportion crossed each boundary) data.  This is the value that will be miniized by the optimization program.
#' @keywords RRW random walk assess fit
#' @export
#' @examples assessRRWfit (data, b=14, s=0.1, loopsPerRWstep = 400)

assessRRWfit <- function(data, b, s=0, nSD=0, db=0, da=0.2, vc = 0, bCols = NULL, sCols = NULL, nSDCols = NULL, dbCols = NULL, daCols = NULL, vcCols = NULL, dataOverlapCol = "overlap", RwSamplesCol = "Q50", dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = 200, minimizeStat = 'BIC', pars.n = NULL, equalizeRTandPhit = FALSE) {

  #make sure minimizeStat is valid
  minimizeOpts <- c("BIC", "AIC", "R_Square")
  if(!(minimizeStat %in% minimizeOpts) ) {
    stop (paste("you set minimizeStat to:", minimizeStat, ", but it must be one of the following:", minimizeOpts, sep=" "))
  }

  #make sure that the input parameters from the grid search are valid
  validParams <- validateParameters(data=data, b=b, s=s, nSD=nSD, db=db, da=da, vc=vc, bCols=bCols, sCols=sCols, nSDCols=nSDCols, dbCols=dbCols, daCols=daCols, vcCols=vcCols, dataOverlapCol=dataOverlapCol)

  #if the parameters from the grid search are valid, see how well they fit the data
  if(validParams == TRUE) {
    #get the keep and merge columns
    #the RW columns are defined by the program. The column used to compare with the data RT is a choice of the user.
    RWkeepColumns <- c("overlap", RwSamplesCol, "pCross", "correct")
    mergeByDataColumns <- c(dataOverlapCol, dataCorrectCol)
    mergeByRWColumns <- c("overlap", "correct")

    #add the parameter effect columns. if they are NULL, they don't appear anyways
    RWkeepColumns <- c(RWkeepColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols)
    mergeByDataColumns <- c(mergeByDataColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols)
    mergeByRWColumns <- c(mergeByRWColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols)

    df.fitted <- getPredictedRRWpoints(data = data, RWkeepColumns = RWkeepColumns, mergeByDataColumns = mergeByDataColumns, mergeByRWColumns = mergeByRWColumns, dataRtCol = dataRtCol, RwSamplesCol = RwSamplesCol, dataOverlapCol = dataOverlapCol, b = b, startValue = s,  noiseSD = nSD, decayBeta = db, decayAsymptote = da, valueChange = vc, bCols = bCols, sCols = sCols, nSDCols = nSDCols, dbCols = dbCols, daCols = daCols, vcCols = vcCols, loops = loopsPerRWstep)

    #get potential minimization variable (1-r2)
    out.rss <- getMinR2RRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], equalizeRTandPhit = equalizeRTandPhit)

    #get potential minimization variable BIC
    if(is.null(pars.n)) {
      pars.n = 1 + length(b) + ifelse(length(s)==1 & s[1]==0, 0, length(s)) + ifelse(length(nSD)==1 & nSD[1]==0, 0, length(nSD)) + ifelse(length(db)[1]==1 & db==0, 0, length(db)) + ifelse(length(da)[1]==1 & da==0.2, 0, length(da)) + ifelse(length(vc)==1 & vc[1]==0, 0, length(vc))
    }
    out.BIC <- getICRRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "BIC")
    out.AIC <- getICRRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "AIC")

  } else {
    #if the parameters are not valid, then assign the minimization variable Infinite value (Inf)
    out.rss <- Inf
    out.BIC <- Inf
    out.AIC <- Inf
  }

  #output the minimization variable
  switch(minimizeStat, BIC = return(out.BIC), AIC = return(out.AIC), R_Square = return(out.rss))
}
