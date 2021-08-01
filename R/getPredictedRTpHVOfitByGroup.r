#' A function to fit the RT and pHVO function to choice data for different groups
#'
#' This function attempts to fit a non-linear, decelerating function to the pHVO data; a linear function to the RT data.
#' @param data A dataframe containing the choice data (often after running through ch.moralsDataPrep()).
#' @param grpCols A vector of strings containing the grouping columns.  If NULL, then the fit will be run without grouping. DEFAULT = NULL.
#' @param RTcol a string that specifies the name of the column in "data" that contains the RT for each trial.
#' @param overlapRoundCol a string that specifies the name of the column in "data" that contains the overlap column.
#' @param minNperOverlap an integer that specifies the minimum number of trials necessary to include an overlap bin in the graph. DEFAULT = 0.
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = TRUE.
#' @keywords  RT pHit fit plots
#' @return a list containing: the input data and best fit for RT and pHVO; RTfit = the lm object for each group; pHitFit = the nls fit object with the fit of the pHit data for each group.
#' @export
#' @examples getPredictedRTpHVOfitForModelList (data=moralsData,simpleModelList, "resdRT", "overlapRound", minNperOverlap = 40, useTwoParameterModel = TRUE)

getPredictedRTpHVOfitByGroup <- function (data, grpCols = NULL,RTcol, overlapRoundCol, minNperOverlap = 0, useTwoParameterModel = TRUE) {

	regFits <- NULL
	RTfit <- list()
	pHVOFit <- list()
	if(!is.null(grpCols)) {
		df.list <- chutils::ch.subsetDFbyGroups(data, grpCols)
		nDFs <- nrow(df.list$dfIndex)
		#for each group
		for(i in 1:nDFs) {
			df.tmp <- df.list[[i]][df.list[[i]][["correct"]] & df.list[[i]][["n"]] > minNperOverlap,]
			fits.out <- getPredictedRTpHVOfit(df.tmp, "rt", "pHit", overlapRoundCol, useTwoParameterModel = useTwoParameterModel)
			pHVOFit[[i]] <- fits.out$pHVOFit
			RTfit[[i]] <- fits.out$RTfit
			regFits <- chutils::ch.rbind(regFits,fits.out$data)
		}
	} else {
		df.tmp <- data[data[["correct"]] & data[["n"]] > minNperOverlap,]
		fits.out <- getPredictedRTpHVOfit(df.tmp, "rt", "pHit", overlapRoundCol, useTwoParameterModel = useTwoParameterModel)
		pHVOFit[[1]] <- fits.out$pHVOFit
		RTfit[[1]] <- fits.out$RTfit
		regFits <- chutils::ch.rbind(regFits,fits.out$data)
	}
  return(list(data = regFits, RTfit = RTfit, pHVOFit = pHVOFit))
}
