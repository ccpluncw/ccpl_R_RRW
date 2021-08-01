#' A function to fit the RT and pHVO function to choice data
#'
#' This function attempts to fit a non-linear, decelerating function to the pHVO data; a linear function to the RT data.
#' @param data A dataframe containing the summarized choice data.
#' @param RTcol a string that specifies the name of the column in "data" that contains the RT for each cell.
#' @param pHVOcol a string that specifies the name of the new column that will contains a p(HVO) for each cell.
#' @param overlapRoundCol a string that specifies the name of the column in "data" that contains the overlap column.
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = FALSE.
#' @keywords  RT pHit fit plots
#' @return a list containing: the input data and best fit for RT and pHVO; RTfit = the lm object; pHitFit = the nls fit object with the fit of the pHVO data.
#' @export
#' @examples getPredictedRTpHVOfit (data=moralsData,"resdRT", "pHit", "overlapRound", useTwoParameterModel = FALSE)

getPredictedRTpHVOfit <- function (data, RTcol, pHVOcol, overlapRoundCol,  useTwoParameterModel = FALSE) {

	pHVOFit <- NULL
	RTfit <- NULL
	if(!is.null(data)) {
		#set the fits to NA
		data$pHVOfit <- NA
		data$RTfit <- NA
		#make sure there are at least 3 x-axis categories
		if(length(data[[overlapRoundCol]]) > 2) {
			pHVOFit <- chutils::ch.pHVOfit(data[[overlapRoundCol]], data[[pHVOcol]], useTwoParameterModel = useTwoParameterModel)
			if(!is.null(pHVOFit)) {
				data$pHVOfit <- fitted(pHVOFit[["nlsObject"]])
			}
			RTfit <- lm(data[[RTcol]]~data[[overlapRoundCol]])
			if(!is.null(RTfit)) {
				data$RTfit <- fitted(RTfit)
			}
		}
	}

	return(list(data = data, pHVOFit = pHVOFit, RTfit = RTfit))
}
