#' A function to fit the RT and pHVO function to choice data for different groups
#'
#' This function attempts to fit a non-linear, decelerating function to the pHVO data; a linear function to the RT data.
#' @param data A dataframe containing the choice data (often after running through ch.moralsDataPrep()).
#' @param grpCols A vector of strings containing the grouping columns. Do not include the column that identifies the whether the stimuli are above or below the reference distribution value here.  If NULL, then the fit will be run without grouping. DEFAULT = NULL.
#' @param RTcol a string that specifies the name of the column in "data" that contains the RT for each trial.
#' @param pHVOcol a string that specifies the name of the new column that will contains a p(HVO) for each cell.
#' @param overlapRoundCol a string that specifies the name of the column in "data" that contains the overlap column.
#' @param refDistCol a string that specifies the name of the column in "data" that identifies those stimuli whose value are above the reference distribution and those that are below the reference distribution. When this variable is included and useTwoParameterModel is set to TRUE, then the alpha parameter of the p(HVO) function will be fit so it is symetric above and below 0.5 for refHVO and refLVO. Also, the intercept parameter of the RT function is allowed to vary for refHVO and refLVO, but the slope is kept constant. DEFAULT = NULL (the refDistCol is ignored)
#' @param correctCol a string that specifies the name of the column that specifies if the participant chose the item with the greatest value distribution (correct) or if they did not (incorrect). This is used to filter the condition so that the analysis is run only on correct trials (p(HVO) of course, contains the influence of both correct and incorrect trials). If this is set to NULL, then no filtering will occur. DEFAULT = NULL
#' @param correctVals a vector of two values that specifies the "correct" value (index 1) and the "incorrect" value (index 2). e.g, c("yes", "no"). DEFAULT = c(TRUE, FALSE)
#' @param nCol a string that specifies the name of the column that specifies the number of trials that composed the condition in each row of the dataframe. This is used to filter the data by minNperOverlap.  If NULL, then there will be no filtering by n. DEFAULT = NULL.
#' @param minNperOverlap an integer that specifies the minimum number of trials necessary to include an overlap bin in the graph. DEFAULT = 0.
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = TRUE.
#' @keywords  RT pHit fit plots
#' @return a list containing: the input data and best fit for RT and pHVO; RTfit = the lm object for each group; pHitFit = the nls fit object with the fit of the pHit data for each group.
#' @export
#' @examples getPredictedRTpHVOfitForModelList (data=moralsData,simpleModelList, "resdRT", "overlapRound", minNperOverlap = 40, useTwoParameterModel = TRUE)

getPredictedRTpHVOfitForRefDist <- function (data, grpCols = NULL,RTcol, pHVOcol, overlapRoundCol, refDistCol = NULL, correctCol = NULL, correctVals = c(TRUE, FALSE), nCol = NULL, minNperOverlap = 0, useTwoParameterModel = TRUE) {

	regFits <- NULL
	RTfit <- list()
	pHVOFit <- list()
	RTfitLVO <- list()
	if(!is.null(grpCols)) {
		df.list <- chutils::ch.subsetDFbyGroups(data, grpCols)
		nDFs <- nrow(df.list$dfIndex)
		#for each group
		for(i in 1:nDFs) {
			df.tmp <- df.list[[i]]
			if(!is.null(nCol)) {
				df.tmp <- df.tmp[df.tmp[[nCol]] > minNperOverlap,]
			}
			fits.out <- getPredictedRTpHVOfit(df.tmp, RTcol, pHVOcol, overlapRoundCol, grpCol = refDistCol, correctCol = correctCol, correctVals = correctVals, useTwoParameterModel = useTwoParameterModel)
			pHVOFit[[i]] <- fits.out$pHVOFit
			RTfit[[i]] <- fits.out$RTfit
			RTfitLVO[[i]] <- fits.out$RTfitLVO
			regFits <- chutils::ch.rbind(regFits,fits.out$data)
		}
	} else {
		df.tmp <- data
		if(!is.null(nCol)) {
			df.tmp <- df.tmp[df.tmp[[nCol]] > minNperOverlap,]
		}
		fits.out <- getPredictedRTpHVOfit(df.tmp, RTcol, pHVOcol, overlapRoundCol, grpCol = refDistCol, correctCol = correctCol, correctVals = correctVals, useTwoParameterModel = useTwoParameterModel)
		pHVOFit[[1]] <- fits.out$pHVOFit
		RTfit[[1]] <- fits.out$RTfit
		RTfitLVO[[1]] <- fits.out$RTfitLVO
		regFits <- chutils::ch.rbind(regFits,fits.out$data)
	}
  return(list(data = regFits, pHVOFit = pHVOFit, RTfit = RTfit, RTfitLVO = RTfitLVO))
}
