#' Function that converts a evaluation criterion shift (in SD units) to an overlap. This is essentially identical to "convertValueShiftToOverlap()"
#'
#' This function converts a evaluation criterion shift (in SD units) to an overlap. To do so, it uses an equation fit to normal distribution changes - which will fit the average distibution shape (if not all individual distributions)
#' @param currentOverlap A number that is the current overlap of the two distributions before the mean shift.
#' @param criterionShiftInSDunits A number that is the evaluation criterion shift shift SD units.
#' @param allowExactZeroToShift A boolean that specifies whether to allow an exact zero overlap to shift. You may want to set this to FALSE if you beleive that an overlap of Zero indicated distributions that are so far appart that any mean shift if unlikely to change that overlap. When it is set to TRUE, the overlap of Zero is set to 0.0001 (allowing a shift to occur). This only applies to "exact zeros." When this function reduces an overlap to Zero, that reduction size is known and will not require a correction. DEFAULT = TRUE.
#' @return A number that represents the new overlap after the mean shift is taken into account
#' @keywords overlap shift mean vc value change
#' @export
#' @examples convertCriterionShiftToOverlap (0.5, 1)

convertCriterionShiftToOverlap <- function (currentOverlap, criterionShiftInSDunits, allowExactZeroToShift = T) {
	overOneCorrection <- FALSE
	#  This is aprt of the overOneCorrection routine
	criterionShiftSign <- sign(criterionShiftInSDunits)

	newOverlap <- currentOverlap

	### only run the analysis if the criterion actually shifts
	if(criterionShiftInSDunits != 0 & (newOverlap <= 0 & criterionShiftSign < 0) == FALSE & (newOverlap >= 2 & criterionShiftSign > 0) == FALSE) {
				#if newOverlap >=0
				if(newOverlap >= 0  & newOverlap <= 2) {
					#if newOverlap is greater than one (it will need a bit of correction).
					if(newOverlap > 1 & newOverlap <= 2) {
						#Here subtract newOverlap from 2 so we get an overlap between 0 and 1
						newOverlap <- 2 - newOverlap
						#set flag that the correction has occured
						overOneCorrection <- TRUE
					} else {
						overOneCorrection <- FALSE
					}
					#first allow newOverlap to change if it == 0
					if(allowExactZeroToShift & newOverlap == 0) newOverlap <- 0.0001
					#this is a formula I calculated from shifting two normal distributions relative to one another
					#calculate criterion shift required for the function to fit. Here, we add or subtract a constant to the criterion.
					criterionPreCalcShift <- (5.240*newOverlap^0.2593 - 5.240)

					#### add correction factor that increases the accuracy of the estimation
					cf <- ifelse(newOverlap <= 0.62, 0.00351 + (-0.05533 - 0.00351)/(1 + exp(3.33846 - 10.21863*newOverlap)), -0.05678 + (0.06126 - -0.05678)/(1 + exp(9.25117 - 9.38254*newOverlap)))

					#add criterion shift to the calculation
					#When standard is HVO, then add a criterion shift
					#when standard is LVO, then subtract a criterion shift
					if(overOneCorrection) cShiftSign <- -1
						else cShiftSign <- 1

					c2 <- (cShiftSign * criterionShiftInSDunits) + criterionPreCalcShift + cf

					#now calculate the new Overlap
					#from a fit to an analysis
					newOverlap <- (1.65103 + (-0.01553 - 1.65103)/(1+exp(0.42817 - 1.27936 * abs(c2))))

					#if there was an over one corection, convert the result back to the original scale
					#if the new criterion (c2) is greater than 0, then the distributions have switched sides
					if(c2 > 0| overOneCorrection) {
						newOverlap <- (1 - newOverlap) + 1
					}
			}
	}
	if(newOverlap < 0) newOverlap <- 0
	if(newOverlap > 2) newOverlap <- 2

return(newOverlap)
}
