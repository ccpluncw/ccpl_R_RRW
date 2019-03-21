#' This function fits the RRW and returns (1-r2).
#'
#' Function that fits the RRW and returns (1-r2), where r2 is the fit of the RRW simulation to the empirical RT and error data
#'
#'
#'
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contain a column specifying a condition that will influence either the startpoint or the value of the trials.
#' @param b A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param s A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  s = 0 is the default and represents an unbiased start point.
#' @param nSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. nSD = 0 is the default and represents no noise being added to each step.
#' @param db A signed number representing the beta value in the Information Accrual Bias (IAB). db = 0 is the default and represents no IAB.
#' @param da A positive number representing the asymptote value in the Information Accrual Bias (IAB). da = 0.2 is the default.
#' @param startBiasEffect A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  The startBiasEffect is a constant that is added to s (the start point) as a function of the values contained in the startBiasEffectCol. The default is NULL because there is no startBiasEffect unless one includes a startBiasEffectCol.
#' @param startBiasEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the position of the start point.  This column should be effect coded, because the values contained in this column will be multiplied by the value of startBiasEffect (and then added to s) to get the start point of the RRW.
#' @param valueChangeEffect A signed number between that indicates the change of overlap that one predicts as a function of the values contained in the valueChangeEffectCol. The default is NULL because there is no valueChangeEffect unless one includes a valueChangeEffectCol.
#' @param valueChangeEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the change of overlap (value change).  This column should be effect coded, because the values contained in this column will be multiplied by the value of valueChangeEffect (and then added to overlap) to get the overlap input in the RRW. This change is assumed to be constant across the entire overlap sequence (from 0-1).
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param RwSamplesCol A string that identifies the name of the column in the RRW simulations that contains the summary of the samples that you want to use as a simulation for RT.  The possible columns are: "Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile).  The default is "Q50" because the median is more robust than the mean.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param dataCorrectCol A string that identifies the name of the column in data that identifies whether the trials were correct (TRUE) or incorrect (FALSE). The default is "correct"
#' @param loopsPerRWstep A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#'
#' @return (1-r_square) for the fit of the model to the RT and pHit (proportion crossed each boundary) data.  This is the value that will be miniized by the optimazation program.
#' @keywords RRW random walk assess fit
#' @export
#' @examples assessRRWfit (data, b=14, s=0.1, loopsPerRWstep = 400)

assessRRWfit <- function(data, b, s=0, nSD=0, db=0, da=0.2, startBiasEffect = NULL, startBiasEffectCol = NULL, valueChangeEffect = NULL, valueChangeEffectCol = NULL, dataOverlapCol = "overlap", RwSamplesCol = "Q50", dataRtCol = "rt", dataPhitCol = "pHit", dataCorrectCol = "correct", loopsPerRWstep = 200) {
  overlapSeq <- unique(data$overlap)

  #make sure that the input parameters from the grid search are valid
  validParams <- validateParameters(data=data, b=b, s=s, nSD=nSD, db=db, da=da, startBiasEffect = startBiasEffect, startBiasEffectCol = startBiasEffectCol, valueChangeEffect = valueChangeEffect, valueChangeEffectCol = valueChangeEffectCol)

  #if the parameters from the grid search are valid, see how well they fit the data
  if(validParams == TRUE) {

    #get the keep and merge columns
    #the RW columns are defined by the program. The column used to compare with the data RT is a choice of the user.
    RWkeepColumns <- c("overlap", RwSamplesCol, "pCross", "correct")
    mergeByDataColumns <- c(dataOverlapCol, dataCorrectCol)
    mergeByRWColumns <- c("overlap", "correct")

    #add the startBiasEffectCol if it's not NULL
    if(!is.null(startBiasEffectCol)) {
      RWkeepColumns <- c(RWkeepColumns, startBiasEffectCol)
      mergeByDataColumns <- c(mergeByDataColumns, startBiasEffectCol)
      mergeByRWColumns <- c(mergeByRWColumns, startBiasEffectCol)
    }
    #add the ValueChangeEffectCol if it's not NULL
    if(!is.null(valueChangeEffectCol)) {
      RWkeepColumns <- c(RWkeepColumns, valueChangeEffectCol)
      mergeByDataColumns <- c(mergeByDataColumns, valueChangeEffectCol)
      mergeByRWColumns <- c(mergeByRWColumns, valueChangeEffectCol)
    }

    df.fitted <- getPredictedRRWpoints(data, b, RWkeepColumns, mergeByDataColumns, mergeByRWColumns, dataRtCol, RwSamplesCol, startValue = s,  noiseSD = nSD, decayAsymptote = da, decayBeta = db, startBiasEffect = startBiasEffect, startBiasEffectCol = startBiasEffectCol, valueChangeEffect = valueChangeEffect, valueChangeEffectCol = valueChangeEffectCol, loops = loopsPerRWstep, progress = F)

    #get the  minimization variable (1-r2) for the pHit
    pCor.rss <- min.r2(df.fitted[df.fitted$correct == TRUE, "pCross"], df.fitted[df.fitted$correct == TRUE, dataPhitCol])

    #try to regress Q50 on rt.  If it works, get (1-r.squared). Otherwise set (1-r.squared) to 1.
    Q50.rss <- min.r2(df.fitted$rtFit,df.fitted[[dataRtCol]])

    #combine minimization variable (1-r2) for pHit and RT
    out.rss <- (Q50.rss + pCor.rss)/2

  } else {

    #if the parameters are not valid, then assign the minimization variable (1-r2) = 1
    out.rss <- 1

  }

  #output the minimization variable (1-r2)
  return(out.rss)
}
