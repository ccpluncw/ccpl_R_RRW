#' This function validates the parameters for RRW.
#'
#' This is an internal function that validates the parameters for RRW.  It is used with assessRRWfit(), which is used with an opimization program.  Validate parameters ensures that the parameters input are valid. There are no default values because all values must be assessed for validity.
#'
#'
#'
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contain a column specifying a condition that will influence either the startpoint or the value of the trials.
#' @param b A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param s A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.
#' @param nSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step.
#' @param db A signed number representing the beta value in the Information Accrual Bias (IAB).
#' @param da A positive number representing the asymptote value in the Information Accrual Bias (IAB).
#' @param startBiasEffect A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  The startBiasEffect is a constant that is added to s (the start point) as a function of the values contained in the startBiasEffectCol.
#' @param startBiasEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the position of the start point.  This column should be effect coded, because the values contained in this column will be multiplied by the value of startBiasEffect (and then added to s) to get the start point of the RRW.
#' @param valueChangeEffect A signed number between that indicates the change of overlap that one predicts as a function of the values contained in the valueChangeEffectCol.
#' @param valueChangeEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the change of overlap (value change).  This column should be effect coded, because the values contained in this column will be multiplied by the value of valueChangeEffect (and then added to overlap) to get the overlap input in the RRW. This change is assumed to be constant across the entire overlap sequence (from 0-1).
#'
#' @return TRUE if the parameters are valid; FALSE if the parameters are invalid
#' 
#' @keywords RRW random walk parameters validate
#' @export
#' @examples validateParameters (data, b=14, s=0.1, nSD = 0, db = .2, da = .2, startBiasEffect = NULL, startBiasEffectCol = NULL, valueChangeEffect = NULL, valueChangeEffectCol = NULL)

validateParameters <- function(data, b, s, nSD, db, da, startBiasEffect, startBiasEffectCol, valueChangeEffect, valueChangeEffectCol) {

  out <- TRUE

  if(is.null(startBiasEffect)) {
    sIn <- s
  } else {
    if(is.null(startBiasEffectCol)) {
      out <- FALSE
      sIn <- s + 1 * startBiasEffect
    } else {
      sbConds <- unique(data[[startBiasEffectCol]])
      sIn <- s + sbConds * startBiasEffect
    }
  }

  if(is.null(valueChangeEffect)) {
  } else {
    if(is.null(valueChangeEffectCol)) {
      out <- FALSE
    } else {
      vcConds <- unique(data[[valueChangeEffectCol]])
      out <- ifelse( max(abs(vcConds * valueChangeEffect)) >1, FALSE, TRUE)
    }
  }

  if( max(abs(sIn)) >=1 ) {
    out <- FALSE
  }

  if(da < 0 | da >= 1) {
    out <- FALSE
  }

  if(b <= 1) {
    out <- FALSE
  }

  return(out)

}
