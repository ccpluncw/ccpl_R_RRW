#' A function to fit the RT and pHVO function to choice data
#'
#' This function attempts to fit a non-linear, decelerating function to the pHVO data; a linear function to the RT data.
#' @param data A dataframe containing the summarized choice data.
#' @param RTcol a string that specifies the name of the column in "data" that contains the RT for each cell.
#' @param pHVOcol a string that specifies the name of the new column that will contains a p(HVO) for each cell.
#' @param overlapRoundCol a string that specifies the name of the column in "data" that contains the overlap column.
#' @param grpCol a string that specifies the name of the column in "data" that contains the grouping variable that identifies those stimuli whose value are above the reference distribution and those that are below the reference distribution. When this variable is included and useTwoParameterModel is set to TRUE, then the alpha parameter of the p(HVO) function will be fit so it is symetric above and below 0.5 for refHVO and refLVO. Also, the intercept parameter of the RT function is allowed to vary for refHVO and refLVO, but the slope is kept constant. DEFAULT = NULL (the grouping variable is ignored)
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = FALSE.
#' @keywords  RT pHit fit plots
#' @return a list containing: the input data and best fit for RT and pHVO; RTfit = the lm object; pHitFit = the nls fit object with the fit of the pHVO data.
#' @export
#' @examples getPredictedRTpHVOfit (data=moralsData,"resdRT", "pHit", "overlapRound", grpCol = "refValue", useTwoParameterModel = FALSE)

getPredictedRTpHVOfitOld <- function (data, RTcol, pHVOcol, overlapRoundCol, grpCol = NULL, useTwoParameterModel = FALSE) {

	pHVOFit <- NULL
	RTfit <- NULL
	if(!is.null(data)) {
		#set the fits to NA
		data$pHVOfit <- NA
		data$RTfit <- NA
		#make sure there are at least 3 x-axis categories
		if(length(data[[overlapRoundCol]]) > 2) {
			if(is.null(grpCol)) {
				pHVOFit <- chutils::ch.pHVOfit(data[[overlapRoundCol]], data[[pHVOcol]], grp = NULL, useTwoParameterModel = useTwoParameterModel)
				RTfit <- chutils::ch.RTfit(data[[overlapRoundCol]], data[[RTcol]], grp = NULL)
			} else {
				pHVOFit <- chutils::ch.pHVOfit(data[[overlapRoundCol]], data[[pHVOcol]], data[[grpCol]], useTwoParameterModel = useTwoParameterModel)
				RTfit <- chutils::ch.RTfit(data[[overlapRoundCol]], data[[RTcol]], data[[grpCol]])
			}
			if(!is.null(pHVOFit)) {
				data$pHVOfit <- fitted(pHVOFit[["nlsObject"]])
			}
			if(!is.null(RTfit)) {
				data$RTfit <- fitted(RTfit[["RTObject"]])
			}
		}
	}

	return(list(data = data, pHVOFit = pHVOFit, RTfit = RTfit))
}
