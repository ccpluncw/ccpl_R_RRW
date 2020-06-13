#' This function runs the RRW simulation "numSimsToAverage" number of times, each time it fits the simulation to empirical data, and then it averages the fits and outputs the average fit.
#'
#' Function that runs the RRW simulation "numSimsToAverage" number of times, each time it fits the simulation to empirical data, and then it averages the fits and outputs the average fit.  It is generally run after the assessRRWfit() is fit with an optimation routine (e.g., smartGridSearch()) to identify the optimal parameters.  The optimation routine should output the parameters that generate the best fit, and those parameters should be input into plotRRWfit().  This takes the same parameters as assessRRWfit().
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
#' @param numSimsToAverage A number specifying how many times the simulation should be run. When it is run multiple times, the function will calculate the average simulated RT and pHit for each row, and the average and sd r2 and BIC for the set of runs. Default is 10.
#' @param sinkFilename A string that identifies the name of file (.txt) in which the fit statistics will be saved. The default is NULL, whereby the fit statistics will not be saved.
#' @param pars.n A number specifying the number of free parameters in the model.  Use this when you are fixing the values of some parameters, but those fixed values are not the default values of the parameters.  Otherwise, the number of free parameters will be automatically caluculated. Default is NULL (automatically calculate).
#' @param equalizeRTandPhit A boolean that specifies whether the influence of the pHit should be equal to that of rt.  Influence is a function of the number of observations.  RT has more observations than pHit because it has both correct RTs and incorrect RTs.  If this is set to TRUE, then the influence of the pHit and RT is equalized in the minimization statistic. If it is set to FALSE, then the the minimazation statistic is calculated as usual. DEFAULT = FALSE.
#''
#' @return A dataframe that contains the "data" plus the fitted values from the model ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the fitted pHit from the model - probability of crossing each boundary); rtFit (the fitted rt values from the model))
#' @keywords RRW random walk plot output fit
#' @export
#' @examples getRRWfit (data, b=14, s=0.1, loopsPerRWstep = 2000, sinkFilename = "outStats.txt")

getMeanRRWfit <- function(data, b, s=0, nSD=0, db=0, da=0.2, vc = 0, bCols = NULL, sCols = NULL, nSDCols = NULL, dbCols = NULL, daCols = NULL, vcCols = NULL, dataOverlapCol = "overlap", RwSamplesCol = "Q50", dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = 2000, numSimsToAverage = 10, sinkFilename = NULL, pars.n = NULL, equalizeRTandPhit = FALSE) {

  #the RW columns are defined by the program. The column used to compare with the data RT is a choice of the user.
  RWkeepColumns <- c("overlap", RwSamplesCol, "pCross", "correct")
  mergeByDataColumns <- c(dataOverlapCol, dataCorrectCol)
  mergeByRWColumns <- c("overlap", "correct")

  #add the parameter effect columns. if they are NULL, they don't appear anyways
  RWkeepColumns <- c(RWkeepColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols)
  mergeByDataColumns <- c(mergeByDataColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols)
  mergeByRWColumns <- c(mergeByRWColumns, bCols, sCols, nSDCols, dbCols, daCols, vcCols)

  if(is.null(pars.n)) {
    pars.n = 1 + length(b) + ifelse(length(s)==1 & s[1]==0, 0, length(s)) + ifelse(length(nSD)==1 & nSD[1]==0, 0, length(nSD)) + ifelse(length(db)[1]==1 & db==0, 0, length(db)) + ifelse(length(da)[1]==1 & da==0.2, 0, length(da)) + ifelse(length(vc)==1 & vc[1]==0, 0, length(vc))
  }

  df.fitted <- NULL
  pCor.rss <- NULL
  Q50.rss <- NULL
  ave.rss <- NULL
  out.rss <- NULL
  out.BIC <- NULL
  out.AIC <- NULL
  simPhitCols <- NULL
  simSMPSCols <- NULL
  simRTCols <- NULL
  for (i in 1:numSimsToAverage) {

    cat(i," of ", numSimsToAverage, "\n")

    df.tmp <- getPredictedRRWpoints(data = data, RWkeepColumns = RWkeepColumns, mergeByDataColumns = mergeByDataColumns, mergeByRWColumns = mergeByRWColumns, dataRtCol = dataRtCol, RwSamplesCol = RwSamplesCol, dataOverlapCol = dataOverlapCol, b = b, startValue = s,  noiseSD = nSD, decayBeta = db, decayAsymptote = da, valueChange = vc, bCols = bCols, sCols = sCols, nSDCols = nSDCols, dbCols = dbCols, daCols = daCols, vcCols = vcCols, loops = loopsPerRWstep)

    #get (1-r2) for the pHit and RT.
     pCor.rss[i] <- getMinR2(df.tmp[df.tmp$correct == TRUE, "pCross"], df.tmp[df.tmp$correct == TRUE, dataPhitCol], na.rm = T)
     Q50.rss[i] <- getMinR2(df.tmp$rtFit,df.tmp[[dataRtCol]], na.rm = T)
    # #combine (1-r2) for pHit and RT
    ave.rss[i] <- (Q50.rss[i] + pCor.rss[i])/2

    out.rss[i] <- getMinR2RRW(df.tmp[df.tmp$correct == TRUE, "pCross"], df.tmp[df.tmp$correct == TRUE, dataPhitCol], df.tmp$rtFit,df.tmp[[dataRtCol]], equalizeRTandPhit = equalizeRTandPhit)

    #calculate BIC
    out.BIC[i] <- getICRRW(df.tmp[df.tmp$correct == TRUE, "pCross"], df.tmp[df.tmp$correct == TRUE, dataPhitCol], df.tmp$rtFit,df.tmp[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "BIC")
    out.AIC[i] <- getICRRW(df.tmp[df.tmp$correct == TRUE, "pCross"], df.tmp[df.tmp$correct == TRUE, dataPhitCol], df.tmp$rtFit,df.tmp[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "AIC")


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

  #rescale so samples fits rts
  Q50.lm <- lm(df.fitted[[dataRtCol]]~df.fitted[[RwSamplesCol]])
  df.fitted$rtFit = df.fitted[[RwSamplesCol]]*coef(Q50.lm)[2]+coef(Q50.lm)[1]

  #caluculate (1-r2) on the average of fits the pHit and RT.
   pCor.rss.final <- getMinR2(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], na.rm=T)
   Q50.rss.final <- getMinR2(df.fitted$rtFit,df.fitted[[dataRtCol]], na.rm=T)
  # #combine (1-r2) for pHit and RT
  ave.rss.final <- (Q50.rss.final + pCor.rss.final)/2

  out.rss.final <- getMinR2RRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], equalizeRTandPhit = equalizeRTandPhit)

  #calculate BIC on the average of fits
  out.BIC.final <- getICRRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "BIC")
  out.AIC.final <- getICRRW(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol], df.fitted$rtFit,df.fitted[[dataRtCol]], pars.n, equalizeRTandPhit = equalizeRTandPhit, ICtype = "AIC")

  #get mean and sd of fit stats of the individual runs
  m.Q50.rss <- mean(Q50.rss)
  sd.Q50.rss <- sd(Q50.rss)
  m.pCor.rss <- mean(pCor.rss)
  sd.pCor.rss <- sd(pCor.rss)
  m.ave.rss <- mean(ave.rss)
  sd.ave.rss <- sd(ave.rss)
  m.out.rss <- mean(out.rss)
  sd.out.rss <- sd(out.rss)
  m.out.BIC <- mean(out.BIC)
  sd.out.BIC <- sd(out.BIC)
  m.out.AIC <- mean(out.AIC)
  sd.out.AIC <- sd(out.AIC)

  if (!is.null(sinkFilename)) {
    sink(sinkFilename)
      cat("\n\n ******** Simulation Epoch Transformation ******** \n\n")
      cat(" RTb = ", coef(Q50.lm)[2], "\n\n")

      cat("\n\n ******** Parameter Values ******** \n\n")
      cat(" RTi (Ter) = ", coef(Q50.lm)[1], "\n\n")

      if(!is.null(bCols)) cat(" Boundary Effect Codes = ", bCols, "\n")
      cat(" Boundary = ", b, "\n\n")

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

      cat("\n\n ******** pHit and RT R Square ******** \n\n")
      cat(" Number of Simulations Run with the parameter values = ", numSimsToAverage, "\n")
      cat("\n mean p(hit) R Square = ", 1 - m.pCor.rss, "\n")
      cat(" sd p(hit) R Square = ", sd.pCor.rss, "\n")
      cat(" Final p(hit) R Square = ", 1 - pCor.rss.final, "\n")
      cat("\n mean RT R Square = ", 1 - m.Q50.rss, "\n")
      cat(" sd RT R Square = ", sd.Q50.rss, "\n")
      cat(" Final RT R Square = ", 1 - Q50.rss.final, "\n")
      cat("\n mean Average R Square = ", 1 - m.ave.rss, "\n")
      cat(" sd AverageR Square = ", sd.ave.rss, "\n")
      cat(" Final Average R Square = ", 1 - ave.rss.final, "\n")

      cat("\n\n ******** Model Fit Statistics ******** \n\n")
      cat(" N Parameters = ", pars.n, "\n")
      cat("\n Average Overall R Square = ", 1-m.out.rss, "\n")
      cat(" SD Overall R Square = ", sd.out.rss, "\n")
      cat(" Final Overall R Square = ", 1 - out.rss.final, "\n")
      cat("\n Average BIC = ", m.out.BIC, "\n")
      cat(" SD BIC = ", sd.out.BIC, "\n")
      cat(" Final BIC = ", out.BIC.final, "\n")
      cat("\n Average AIC = ", m.out.AIC, "\n")
      cat(" SD AIC = ", sd.out.AIC, "\n")
      cat(" Final AIC = ", out.AIC.final, "\n")
      cat("\n Equalize pHit and RT = ", equalizeRTandPhit, "\n\n")
    sink(NULL)
  }

  return(df.fitted)
}
