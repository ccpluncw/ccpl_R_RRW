#' Function that converts a mean shift of one distribution (in SD units) to an overlap. This is statistically the same as "convertCriterionShiftToOverlap()"
#'
#' This function converts a value shift of the LVO distribution (in SD units) to an overlap. To do so, it assumes a two normal distributions with equal variance - which will fit the average distibution shape (if not all individual distributions)
#' @param currentOverlap A number that is the current overlap of the two distributions before the mean shift.
#' @param valueShiftInSDunits A number that is the mean shift of one distribution shift SD units.
#' @param allowExactZeroToShift A boolean that specifies whether to allow an exact zero overlap to shift. You may want to set this to FALSE if you beleive that an overlap of Zero indicated distributions that are so far appart that any mean shift if unlikely to change that overlap. When it is set to TRUE, the overlap of Zero is set to 0.001 (allowing a shift to occur). This only applies to "exact zeros." When this function reduces an overlap to Zero, that reduction size is known and will not require a correction. DEFAULT = TRUE.
#' @return A number that represents the new overlap after the distribution shift is taken into account
#' @keywords overlap shift mean vc value change
#' @export
#' @examples convertDistributionShiftToOverlap (0.5, 1)

convertDistributionShiftToOverlap <- function (currentOverlap, shiftInSDunits, allowExactZeroToShift = T) {
	#first allow currentOverlap to change if it == 0
	if(allowExactZeroToShift & currentOverlap == 0) currentOverlap <- 0.001
	#This corrects for overlaps greater than 1
	overOneCorrection <- FALSE
	
	#get the p associated with current Overlap
	currentP <- chValues::ch.overlapToP(currentOverlap)
	#convert currentP to the difference in distributions in SD units
	currentDiffInSDUnits <- -1*qnorm(currentP, 0,sqrt(2))

	#if the shiftInSDunits will make overlap greater than 1, then do the overOneCorrection
	if(currentDiffInSDUnits + shiftInSDunits  >0 ) overOneCorrection <- TRUE

	#get p associated with the shift in SD units
	pX1GreaterThanX2 <- 1-pnorm(currentDiffInSDUnits, -1*shiftInSDunits, sqrt(2))
	#convert p back to overlap
	newOverlap <- chValues::ch.pToOverlap(pX1GreaterThanX2)

	#if necessary, implement the Over One Correction
	if(overOneCorrection) newOverlap <- (1 - newOverlap) + 1
	
	return(newOverlap)
}

