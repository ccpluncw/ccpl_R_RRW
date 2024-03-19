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
	## Identify the number of iterations the correction must have
	loops <- meanShiftInSDunits/0.02

	#number of iterations need to be positive - 0.02 SD shift is the max accuracy
	units <- abs(trunc(loops))

	## starting values
	newOverlap <- currentOverlap
	diffOverlap <- 0
	# this is a booliean that specifies whether the newOverlap > 1 and needed a correction
	overOneCorrection <- FALSE
	#  This is aprt of the overOneCorrection routine
	meanShiftSign <- sign(meanShiftInSDunits)
	#successively calculate the change in overlap as a function of mean shift in SD units
	while((units) > 0) {
		newOverlapNminus1 <- newOverlap
		#set this iteration of localMeanShiftSign back to the original meanShiftSign
		localMeanShiftSign <- meanShiftSign

		#end loop: if newOverlap == 0 and you are still inreasing the value of the HVO (you can't spread apart when overlap is 0)
		if(newOverlap <= 0 & meanShiftSign > 0) {
			#make new overlap == 0 (just in case it is < 0) and break
			newOverlap <- 0
			break
		}
		#end loop: if newOverlap == 2 then there is 0 overlap, but the LVO became the HVO with the shift. If shift continues to reduce the value of the HVO, then it will actually spread distributions apart. Break here because you can't spread apart when overlap is 0.
		if(newOverlap >= 2 & meanShiftSign < 0) {
			#make new overlap == 2 (just in case it is > 2) and break
			newOverlap <- 2
			break
		}
			#if newOverlap >=0
			if(newOverlap >= 0  & newOverlap <= 2) {
				#if newOverlap is greater than one (it will need a bit of correction).
				if(newOverlap > 1 & newOverlap <= 2) {
					#Here subtract newOverlap from 2 so we get an overlap between 0 and 1
					newOverlap <- 2 - newOverlap
					#set flag that the correction has occured
					overOneCorrection <- TRUE
					#change the sign of this iterations mean shift
					localMeanShiftSign <- localMeanShiftSign*-1
				} else {
					overOneCorrection <- FALSE
				}
				#first allow newOverlap to change if it == 0
				if(allowExactZeroToShift & newOverlap == 0) newOverlap <- 0.0001
				#this is a formula I calculated from shifting two normal distributions relative to one another
				if(localMeanShiftSign < 0) {
					newOverlap <- (1.012 * newOverlap^0.993)
				} else {
					newOverlap <- (0.9884 * newOverlap^1.0091)
				}
				if(overOneCorrection) {
					#if there was an over one corection, convert the result back to the original scale
					newOverlap <- (1 - newOverlap) + 1
				}
			} else {
				newOverlap <- 0
			}
			units <- units - 1
	}
	return(newOverlap)
}
