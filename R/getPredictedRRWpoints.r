#' This function simulates the RRW for a single set of parameters across all values of overlap.
#'
#' This function simulates the RRW for a single set of parameters across all values of overlap and then attaches this simulation to the empirical RT and error data.  This simulation's fit will be assessed by the optimization program (using assessRRWfit()) and then the parameters will be adjested.
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param RWkeepColumns A vector of strings that specify the columns of the RRW simulated data that should be kept when the simulated dataframe is merged with the empirical dataframe.
#' @param mergeByDataColumns A vector of strings that specify the columns of the empirical data that should be used in the "by.x" option of the merge() function when the simulated dataframe is merged with the empirical dataframe.
#' @param mergeByRWColumns A vector of strings that specify the columns of the RRW simulated data that should be used in the "by.y" option of the merge() function when the simulated dataframe is merged with the empirical dataframe.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param RwSamplesCol A string that identifies the name of the column in the RRW simulations that contains the summary of the samples that you want to use as a simulation for RT.  The possible columns are: "Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile).  The default is "Q50" because the median is more robust than the mean.
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row.
#' @param b A vector of number(s) specifying the boundary distance from a 0 startpoint. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bCols" argument. b is specific to the RRW simulation and has no default value
#' @param bcs A vector of positive values that specify the likelihood that the boundary will be reduced because there it is taking too many samples before reaching a threshold. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of b. All effect/dummy codes must >=0. The column names must be specified in the "bcsCols" argument. Essentially, this models a decision-maker reducing the amount of information that they need before a response because it is taking to long.  The boundary reduces when there is little change in the likelihood of a response after N number of steps.  The routine makes this determination as a function of the size of the unaltered boundary (boundaryChangeSensitivity * boundary). Small values (e.g., 0.1), result in an impatient person (relatively early chnage in the boundary values).  Large values (e.g., 0.9), result in a patient person (a boundary that is relatively resistent to change.  A value of 0 is an exception: it will result in no boundary reduction under any circumstances. Impatience increases the likelihood of an error, so it increases the average time of error responses.  Default is 0.25.
#' @param Ter A vector of positive values that adjusts the time of encoding/response (non-decision time). If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter. All effect/dummy codes must >=0. The column names must be specified in the "TerCols" argument. When Ter = 0, a single, best fitting Ter will be calculated. Use this parameter when you have several conditions and you hypothesize that some conditions should have a larger Ter than others.  For example, this may model longer encoding process when different numbers of items are on the screen for different conditions. The hypothetical fastest condition should have Ter = 0, which simply indicates that the RRW program will find the best fit Ter. The slower conditions will have larger Ter values and will be fit ralative to the fastest condition. The change in samples resulting from the Ter effect is a function of the size of the boundary (Ter * boundary). Default = 0, This is still experimental.
#' @param startValue A vector of signed number(s) between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.   If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of startValue. The column names must be specified in the "sCols" argument. startValue = 0 is the default and represents an unbiased start point.
#' @param noiseSD A vector of positive number(s) representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step.  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of noiseSD. The column names must be specified in the "nSDCols" argument. noiseSD = 0 is the default and represents no noise being added to each step.
#' @param decayBeta A vector of signed number(s) representing the beta value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of decayBeta. The column names must be specified in the "dbCols" argument. decayBeta = 0 is the default and represents no IAB.
#' @param decayAsymptote A vector of positive number(s) representing the asymptote value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of decayAsymptote. The column names must be specified in the "daCols" argument. decayAsymptote = 0.2 is the default.
#' @param valueChange A vector of signed number(s) that indicates the change of overlap that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of valueChange. The column names must be specified in the "vcCols" argument. The default is 0 because there is no valueChangeEffect.
#' @param evaluationCriterion A vector of signed number(s) that indicates the change of the evaluation criterion (the position of the criterion in SDT) that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of ec. The column names must be specified in the "ecCols" argument. The default is 0 because this is an unbiased evaluation criterion.
#' @param bCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of boundary (b).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of b.
#' @param bcsCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the boundaryChangeSensitivity (bcs).  These columns must be effect coded, because the values contained in this column will be multiplied by the value of bcs. All effect/dummy codes must be >=0.
#' @param TerCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the non-decision time (Ter).  These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter.  All effect/dummy codes must be >=0.
#' @param sCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of startValue (startValue).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of startValue.
#' @param nSDCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of noiseSD (noiseSD).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of noiseSD.
#' @param dbCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decayBeta (decayBeta).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of decayBeta.
#' @param daCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decay asymptote (decayAsymptote).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of decayAsymptote.
#' @param vcCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of value change (valueChange).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of overlap.
#' @param ecCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of the evaluation criterion (ec).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of overlap.
#' @param loops A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @return A dataframe that contains the "data" plus the fitted values from the model ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the fitted pHit from the model - probability of crossing each boundary); rtFit (the fitted rt values from the model))
#' @keywords RRW random walk simulation get predicted
#' @export
#' @examples getPredictedRRWpoints (data, RWkeepColumns = c("overlap", "Q50", "pCross", "correct"), mergeByDataColumns = c("overlap", "correct"), dataRtCol = "rt", RwSamplesCol = "Q50", dataOverlapCol= "overlap", b=14, s=0.1, loops = 400)

getPredictedRRWpoints <- function (data, RWkeepColumns, mergeByDataColumns, mergeByRWColumns, dataRtCol, RwSamplesCol, dataOverlapCol, b, boundaryChangeSensitivity = 0.25, Ter = 0, startValue = 0,  noiseSD = 0, decayBeta = 0.0, decayAsymptote = 0.2, valueChange = 0, evaluationCriterion = 0, bCols = NULL, bcsCols = NULL, TerCols = NULL, sCols = NULL, nSDCols = NULL, dbCols = NULL, daCols = NULL, vcCols = NULL, ecCols = NULL, loops = 200) {

    #get the unique rows to run the simulation on.  First, you need the overlap column
    dataRunCols <- c(dataOverlapCol, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols, ecCols)
    effectCols <- c(bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols, ecCols)

    #extract the unique information
    data.tmp <- unique(data[,dataRunCols, drop=F])

    data.n <- nrow(data.tmp)
    df.out <- NULL

    #Run the simulation on each unique row
    for(i in 1:data.n) {

      #get the parameter values summed over the columns
#      ovIn <- data.tmp[i,dataOverlapCol] + getRowParameterValue(data.tmp[i,], vcCols, valueChange)

      ######################## START ############################
      ##### When I turn on the new overlap shift calculations, simply uncomment the next two lines, and comment the previous line.
      ##### This code gets mean shift which is a free parameter
      vcMS <- getRowParameterValue(data.tmp[i,], vcCols, valueChange)
      ##### This code converts the mean shift into an overlap
        #this is an old version - never implemented
        #ovIn <- convertValueShiftToOverlap(data.tmp[i,dataOverlapCol], vcMS)
      ovIn <- convertDistributionShiftToOverlap(data.tmp[i,dataOverlapCol], vcMS)
      ##### This code converts the evaluation criterion shift into an overlap. The input is ovIn,
      ##### just in case there is both a value change effect and an evaluation criterion change
      ecMS <- getRowParameterValue(data.tmp[i,], ecCols, evaluationCriterion)
        #this is an old version - never implemented
        #ovIn <- convertCriterionShiftToOverlap(ovIn, ecMS)
      ovIn <- convertDistributionShiftToOverlap(ovIn, ecMS)
      ######################## DONE ############################

      bIn <- getRowParameterValue(data.tmp[i,], bCols, b)
      bcsIn <- getRowParameterValue(data.tmp[i,], bcsCols, boundaryChangeSensitivity)
      TerIn <- getRowParameterValue(data.tmp[i,], TerCols, Ter)
      sIn <- getRowParameterValue(data.tmp[i,], sCols, startValue)
      nSDIn <- getRowParameterValue(data.tmp[i,], nSDCols, noiseSD)
      daIn <- getRowParameterValue(data.tmp[i,], daCols, decayAsymptote)
      dbIn <- getRowParameterValue(data.tmp[i,], dbCols, decayBeta)

      #run the simulation
      df.tmp <- getMomentsOfRRW(overlap = ovIn, b = bIn, boundaryChangeSensitivity = bcsIn, Ter = TerIn, startValue = sIn, noiseSD = nSDIn, decayAsymptote = daIn, decayBeta = dbIn, loops=loops)

      #add the effect codes onto the dataframe containing the simulated data.  You will need that for merging with the actual data
      df.tmp <- cbind(df.tmp, data.tmp[i,effectCols, drop=F])
      df.tmp$overlap <- data.tmp[i,dataOverlapCol]

      #append each simulated row to the output dataframe
      df.out <- chutils::ch.rbind(df.out, df.tmp)
    }

    #combine the simulated data with the empirical data
    df.out1 <- df.out[,RWkeepColumns]

    df.out2 <- merge(df.out1, data, by.x = mergeByRWColumns, by.y=mergeByDataColumns)

    #scale the RRW samples to fit the empirical RT data
    fit.lm <- NULL
    tryCatch ({
      fit.lm <- lm(df.out2[[dataRtCol]]~df.out2[[RwSamplesCol]])
    }, error = function(e) {})
    if(!is.null(fit.lm)) {
      df.out2$rtFit = df.out2[[RwSamplesCol]]*coef(fit.lm)[2]+coef(fit.lm)[1]
    } else {
      df.out2$rtFit = NA
    }

    return(df.out2)

}
