#' Function that converts a mean shift of one distribution (in SD units) to an overlap
#'
#' This function converts a mean shift of one distribution (in SD units) to an overlap. To do so, it uses an equation fit to normal distribution changes - which will fit the average distibution shape (if not all individual distributions)
#' @param currentOverlap A number that is the current overlap of the two distributions before the mean shift. 
#' @param meanShiftInSDunits A number that is the mean shift of one distribution shift SD units.
#' @param allowExactZeroToShift A boolean that specifies whether to allow an exact zero overlap to shift. You may want to set this to FALSE if you beleive that an overlap of Zero indicated distributions that are so far appart that any mean shift if unlikely to change that overlap. When it is set to TRUE, the overlap of Zero is set to 0.0001 (allowing a shift to occur). This only applies to "exact zeros." When this function reduces an overlap to Zero, that reduction size is known and will not require a correction. DEFAULT = TRUE.
#' @return A number that represents the new overlap after the mean shift is taken into account
#' @keywords overlap shift mean vc value change
#' @export
#' @examples convertMeanShiftToOverlap (0.5, 1)

convertMeanShiftToOverlap <- function (currentOverlap, meanShiftInSDunits, allowExactZeroToShift = T) {
	loops <- meanShiftInSDunits/0.05
	dec <- loops%%1
	units <- trunc(loops)
	newOverlap <- currentOverlap
	#successively calculate the change in overlap as a function of mean shift in SD units
	for(i in 1:units) {
		if(newOverlap >= 0) {
			if(allowExactZeroToShift & newOverlap == 0) newOverlap <- 0.0001
				#this is a formula I calculated from shifting two normal distributions relative to one another
			diffOverlap <- (-0.031 * newOverlap^0.61)*sign(meanShiftInSDunits)
			newOverlap <- newOverlap+diffOverlap
		} else {
			newOverlap <- 0
		}
	}
	####### impute decimal version of the change from last difference estimate
	lastDiff <- diffOverlap*dec
	newOverlap <- newOverlap+lastDiff
	return(newOverlap)
}