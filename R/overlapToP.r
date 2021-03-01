#' This function converts Overlap to the p(D1 > D2) .
#'
#' Assume two overlapping distributions, D1 and D2. Randomly select a value from D1 and D2.  p(D1 > D2) is the probability that you select a value from D1 that is greater than D2.  Overlap = 1 - [abs (p(D1 > D2) - 0.5) / 0.5].  This function converts Overlap back to p(D1 > D2).
#' @param overlap A number that equals 1 - [abs (p(D1 > D2) - 0.5) / 0.5].
#' @return p(D1 > D2)
#' @keywords overlap meanP p(D1 > D2)
#' @export
#' @examples overlapToP (0.5)

overlapToP <- function (overlap) {

  #here, as overlap moves greater than 1, that means that D2 starts to become the greater distribution.
  #p gives you the probability of choosing D1.  This vacillates between 0-1 as the overlap changes.

  #first get the modulo of overlap by 4 (this only has an influence if overlap is greater than 4! it should be between 0 and 1)
  ovMod4 <- overlap %% 4
  #then correct the overlap to vacillate properly beween 1 and 0. Again, this only has an influence if overlap is greater than abs(2).
  overlap <- ifelse(ovMod4 > 2, abs(4 - ovMod4), ovMod4)

  p <- (-0.5*(overlap - 1)) + 0.5

  # if(p > 1) {
  #   p <- 1
  # }
  # if(p < 0) {
  #   p <- 0
  # }

  return(p)
}
