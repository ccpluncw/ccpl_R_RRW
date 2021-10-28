#' A function to fit the RT and pHVO function to choice data
#'
#' This function attempts to fit a non-linear, decelerating function to the pHVO data; a linear function to the RT data.
#' @param data A dataframe containing the summarized choice data.
#' @param RTcol a string that specifies the name of the column in "data" that contains the RT for each cell.
#' @param pHVOcol a string that specifies the name of the new column that will contains a p(HVO) for each cell.
#' @param overlapRoundCol a string that specifies the name of the column in "data" that contains the overlap column.
#' @param grpCol a string that specifies the name of the column in "data" that contains the grouping variable that identifies those stimuli whose value are above the reference distribution and those that are below the reference distribution. When this variable is included and useTwoParameterModel is set to TRUE, then the alpha parameter of the p(HVO) function will be fit so it is symetric above and below 0.5 for refHVO and refLVO. Also, the intercept parameter of the RT function is allowed to vary for refHVO and refLVO, but the slope is kept constant. DEFAULT = NULL (the grouping variable is ignored)
#' @param correctCol a string that specifies the name of the column that specifies if the participant chose the item with the greatest value distribution (correct) or if they did not (incorrect). If this is set to NULL, then then a single RT function will be fit to the collapsed correct and incorrect trials. DEFAULT = NULL
#' @param correctVals a vector of two values that specifies the "correct" value (index 1) and the "incorrect" value (index 2). e.g, c("yes", "no"). DEFAULT = c(TRUE, FALSE)
#' @param useTwoParameterModel A boolean that specifies whether to use a two parameter p(HOV) model.  If this is set to TRUE, then this function will fit a p(HVO) model whereby the rightmost point (overlap = 1.0) is not fixed at p(HVO) = 0.5. DEFAULT = FALSE.
#' @keywords  RT pHit fit plots
#' @return a list containing: the input data and best fit for RT and pHVO; RTfit = the lm object; pHitFit = the nls fit object with the fit of the pHVO data.
#' @export
#' @examples getPredictedRTpHVOfit (data=moralsData,"resdRT", "pHit", "overlapRound", grpCol = "refValue", useTwoParameterModel = FALSE)

getPredictedRTpHVOfit <- function (data, RTcol, pHVOcol, overlapRoundCol, grpCol = NULL, correctCol = NULL, correctVals = c(TRUE, FALSE), useTwoParameterModel = FALSE) {

	pHVOFit <- NULL
	RTfit <- NULL
	RTfitLVO <- NULL
	if(!is.null(data)) {
		if(!is.null(correctCol)) {
			data.T <- data[data[[correctCol]] == correctVals[1],]
			data.F <- data[data[[correctCol]] == correctVals[2],]
			data.F$RTfit <- NA
			#must add this so rbind works (columns are the same)
			data.F$pHVOfit <- NA
		} else {
			data.T <- data
			data.F <- NULL
		}
		#set the fits to NA
		data.T$pHVOfit <- NA
		data.T$RTfit <- NA
		#make sure there are at least 3 x-axis categories
		if(length(data.T[[overlapRoundCol]]) > 2) {
			if(is.null(grpCol)) {
				pHVOFit <- chutils::ch.pHVOfit(data.T[[overlapRoundCol]], data.T[[pHVOcol]], grp = NULL, useTwoParameterModel = useTwoParameterModel)
				RTfit <- chutils::ch.RTfit(data.T[[overlapRoundCol]], data.T[[RTcol]], grp = NULL)
				if(!is.null(data.F)) {
					RTfitLVO <- chutils::ch.RTfit(data.F[[overlapRoundCol]], data.F[[RTcol]], grp = NULL)
				}
			} else {
				pHVOFit <- chutils::ch.pHVOfit(data.T[[overlapRoundCol]], data.T[[pHVOcol]], data.T[[grpCol]], useTwoParameterModel = useTwoParameterModel)
				RTfit <- chutils::ch.RTfit(data.T[[overlapRoundCol]], data.T[[RTcol]], data.T[[grpCol]])
				if(!is.null(data.F)) {
					RTfitLVO <- chutils::ch.RTfit(data.F[[overlapRoundCol]], data.F[[RTcol]], data.F[[grpCol]])
				}
			}
			if(!is.null(pHVOFit)) {
				tryCatch ({
					data.T$pHVOfit <- predict(pHVOFit[["nlsObject"]], newdata = data.T)
					}, error = function(e) {
						print(e)
						print(data.T)
				})
			}
			if(!is.null(RTfit)) {
				tryCatch ({
					data.T$RTfit <- predict(RTfit[["RTObject"]], newdata = data.T)
					}, error = function(e) {
						print(e)
						print(data.T)
				})
			}
			if(!is.null(RTfitLVO)) {
				tryCatch ({
					data.F$RTfit <- predict(RTfitLVO[["RTObject"]], newdata = data.F)
					}, error = function(e) {
						print(e)
						print(data.F)
				})
			}
		}
	}

	if(!is.null(data.F)) {
		data <- rbind(data.T, data.F)
	} else {
		data <- data.T
	}
	return(list(data = data, pHVOFit = pHVOFit, RTfit = RTfit, RTfitLVO = RTfitLVO))
}
