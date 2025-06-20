#' This function fits the RRW and returns (1-r2).
#'
#' Function that fits the RRW and returns (1-r2), where r2 is the fit of the RRW simulation to the empirical RT and error data
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param b A vector of number(s) specifying the boundary distance from a 0 startpoint. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bCols" argument. b is specific to the RRW simulation and has no default value
#' @param bcs A vector of positive values that specify the likelihood that the boundary will be reduced because there it is taking too many samples before reaching a threshold. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bcsCols" argument. Essentially, this models a decision-maker reducing the amount of information that they need before a response because it is taking to long.  The boundary reduces when there is little change in the likelihood of a response after N number of steps.  The routine makes this determination as a function of the size of the unaltered boundary (boundaryChangeSensitivity * boundary). Small values (e.g., 0.1), result in an impatient person (relatively early chnage in the boundary values).  Large values (e.g., 0.9), result in a patient person (a boundary that is relatively resistent to change.  A value of 0 is an exception: it will result in no boundary reduction under any circumstances. Impatience increases the likelihood of an error, so it increases the average time of error responses.  Default is 0.25.
#' @param Ter A vector of positive values that adjusts the time of encoding/response (non-decision time). If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter. All effect/dummy codes must >=0. The column names must be specified in the "TerCols" argument. When Ter = 0, a single, best fitting Ter will be calculated. Use this parameter when you have several conditions and you hypothesize that some conditions should have a larger Ter than others.  For example, this may model longer encoding process when different numbers of items are on the screen for different conditions. The hypothetical fastest condition should have Ter = 0, which simply indicates that the RRW program will find the best fit Ter. The slower conditions will have larger Ter values and will be fit ralative to the fastest condition. The change in samples resulting from the Ter effect is a function of the size of the boundary (Ter * boundary). Default = 0, This is still experimental.
#' @param s A vector of signed number(s) between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.   If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of s. The column names must be specified in the "sCols" argument. s = 0 is the default and represents an unbiased start point.
#' @param nSD A vector of positive number(s) representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step.  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD. The column names must be specified in the "nSDCols" argument. nSD = 0 is the default and represents no noise being added to each step.
#' @param db A vector of signed number(s) representing the beta value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of db. The column names must be specified in the "dbCols" argument. db = 0 is the default and represents no IAB.
#' @param da A vector of positive number(s) representing the asymptote value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of da. The column names must be specified in the "daCols" argument. da = 0.2 is the default.
#' @param vc A vector of signed number(s) that indicates the change of overlap that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of vc. The column names must be specified in the "vcCols" argument. The default is 0 because there is no valueChangeEffect.
#' @param ec A vector of signed number(s) that indicates the change of the evaluation criterion (the position of the criterion in SDT) that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of ec. The column names must be specified in the "ecCols" argument. The default is 0 because this is an unbiased evaluation criterion.
#' @param bCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of boundary (b).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of b.
#' @param bcsCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the boundaryChangeSensitivity (bcs).  These columns must be effect coded, because the values contained in this column will be multiplied by the value of bcs. All effect/dummy codes must be >=0.
#' @param TerCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the non-decision time (Ter).  These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter.  All effect/dummy codes must be >=0.
#' @param sCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of startValue (s).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of s.
#' @param nSDCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of noiseSD (nSD).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD.
#' @param dbCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decayBeta (db).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of db.
#' @param daCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decay asymptote (da).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of da.
#' @param vcCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of value change (vc).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of overlap.
#' @param ecCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of the evaluation criterion (ec).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of overlap.
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param RwSamplesCol A string that identifies the name of the column in the RRW simulations that contains the summary of the samples that you want to use as a simulation for RT.  The possible columns are: "Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile).  The default is "Q50" because the median is more robust than the mean.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param dataCorrectCol A string that identifies the name of the column in data that identifies whether the trials were correct (TRUE) or incorrect (FALSE). The default is "correct"
#' @param loopsPerRWstep A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @param minimizeStat A string that specifies which statistic to minimize when optimizing the model fit.  The options are: "BIC" , "AIC" , or "R_Square". Default is "BIC".
#' @param pars.n The number of free parameters. When NULL, the program will attempt to calculate the number of free parameters from the input. Default = NULL.
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.
#' @param numSimsToAverage A number specifying how many times the simulation should be run. When it is run multiple times, the function will calculate the average simulated RT and pHit for each row, and the average and sd r2 and BIC for the set of runs. Default is 10.
#''
#' @return The minimization statistic for the fit of the model to the RT and pHit (proportion crossed each boundary) data.  This is the value that will be miniized by the optimization program.
#' @keywords RRW random walk assess fit
#' @export
#' @examples assessRRWfit (data, b=14, s=0.1, loopsPerRWstep = 400)

assessAveRRWfit <- function(data, b, bcs = 0.25, Ter = 0, s=0, nSD=0, db=0, da=0.2, vc = 0, ec = 0, bCols = NULL, bcsCols = NULL, TerCols = NULL, sCols = NULL, nSDCols = NULL, dbCols = NULL, daCols = NULL, vcCols = NULL, ecCols = NULL, dataOverlapCol = "overlap", RwSamplesCol = "Q50", dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = 200, minimizeStat = 'BIC', pars.n = NULL, equalizeRTandPhit = FALSE, numSimsToAverage = 10) {

  #make sure minimizeStat is valid
  minimizeOpts <- c("BIC", "AIC", "R_Square")
  if(!(minimizeStat %in% minimizeOpts) ) {
    stop (paste("you set minimizeStat to:", minimizeStat, ", but it must be one of the following:", minimizeOpts, sep=" "))
  }

  #make sure that the input parameters from the grid search are valid
  validParams <- validateParameters(data=data, b=b, bcs = bcs, Ter=Ter, s=s, nSD=nSD, db=db, da=da, vc=vc, ec=ec, bCols=bCols, bcsCols=bcsCols, TerCols = NULL, sCols=sCols, nSDCols=nSDCols, dbCols=dbCols, daCols=daCols, vcCols=vcCols, ecCols=ecCols, dataOverlapCol=dataOverlapCol)

  if(is.null(pars.n)) {
    pars.n = 1 + length(b) + ifelse(length(bcs)==1 & bcs[1]==0.25, 0, length(bcs)) + ifelse(length(Ter)==1 & Ter[1]==0, 0, length(Ter)) + ifelse(length(ec)==1 & ec[1]==0, 0, length(ec)) + ifelse(length(s)==1 & s[1]==0, 0, length(s)) + ifelse(length(nSD)==1 & nSD[1]==0, 0, length(nSD)) + ifelse(length(db)[1]==1 & db==0, 0, length(db)) + ifelse(length(da)[1]==1 & da==0.2, 0, length(da)) + ifelse(length(vc)==1 & vc[1]==0, 0, length(vc))
  }

  #if the parameters from the grid search are valid, see how well they fit the data
  if(validParams == TRUE) {
    #get the keep and merge columns
    #the RW columns are defined by the program. The column used to compare with the data RT is a choice of the user.
    RWkeepColumns <- c("overlap", RwSamplesCol, "pCross", "correct")
    mergeByDataColumns <- c(dataOverlapCol, dataCorrectCol)
    mergeByRWColumns <- c("overlap", "correct")

    #add the parameter effect columns. if they are NULL, they don't appear anyways
    RWkeepColumns <- c(RWkeepColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols, ecCols)
    mergeByDataColumns <- c(mergeByDataColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols, ecCols)
    mergeByRWColumns <- c(mergeByRWColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols, ecCols)

    df.fitted <- NULL
    simPhitCols <- NULL
    simSMPSCols <- NULL
    simRTCols <- NULL

    for(i in 1:numSimsToAverage) {
      df.tmp <- getPredictedRRWpoints(data = data, RWkeepColumns = RWkeepColumns, mergeByDataColumns = mergeByDataColumns, mergeByRWColumns = mergeByRWColumns, dataRtCol = dataRtCol, RwSamplesCol = RwSamplesCol, dataOverlapCol = dataOverlapCol, b = b, boundaryChangeSensitivity = bcs,  Ter = Ter, startValue = s,  noiseSD = nSD, decayBeta = db, decayAsymptote = da, valueChange = vc, evaluationCriterion = ec, bCols = bCols, bcsCols = bcsCols, TerCols = TerCols, sCols = sCols, nSDCols = nSDCols, dbCols = dbCols, daCols = daCols, vcCols = vcCols, ecCols = ecCols, loops = loopsPerRWstep)

      #rename the predicted columns to include the run number at the end
      simPhitCols <- c(simPhitCols, paste("pCross",i,sep=""))
      simSMPSCols <- c(simSMPSCols, paste(RwSamplesCol,i,sep=""))
      simRTCols <- c(simRTCols, paste("rtFit",i,sep=""))

      #rename columns
      colnames(df.tmp)[colnames(df.tmp)=="pCross"] <- paste("pCross",i,sep="")
      colnames(df.tmp)[colnames(df.tmp)==RwSamplesCol] <- paste(RwSamplesCol,i,sep="")
      colnames(df.tmp)[colnames(df.tmp)=="rtFit"] <- paste("rtFit",i,sep="")

      #merge the simulated data frame predictions, with the others
      if(is.null(df.fitted)) {
        df.fitted <- df.tmp
      } else {
        df.fitted <- cbind(df.fitted, df.tmp[,c(simPhitCols[i],simSMPSCols[i],simRTCols[i])])
      }
    }

    #get average of simulations
    df.fitted$pCross <- rowMeans(df.fitted[,simPhitCols], na.rm = TRUE)
    df.fitted[[RwSamplesCol]] <- rowMeans(df.fitted[,simSMPSCols], na.rm = TRUE)

    ## remove unstable estimates: Any time less than half the simulations produced a value
    df.Q.NAs <- which(rowMeans(!is.na(df.fitted[,simSMPSCols])) < 0.5)
    df.pc.NAs <- which(rowMeans(!is.na(df.fitted[,simPhitCols])) < 0.5)
    if(length(df.pc.NAs) > 0) {
      df.fitted[df.pc.NAs, "pCross"] <- NA
    }
    if(length(df.Q.NAs) > 0) {
      df.fitted[df.Q.NAs, RwSamplesCol] <- NA
    }

    #rescale so samples fits rts
    Q50.lm <- lm(df.fitted[[dataRtCol]]~df.fitted[[RwSamplesCol]])
    df.fitted$rtFit = df.fitted[[RwSamplesCol]]*coef(Q50.lm)[2]+coef(Q50.lm)[1]

    out.rss <- getMinR2RRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], equalizeRTandPhit = equalizeRTandPhit)
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
