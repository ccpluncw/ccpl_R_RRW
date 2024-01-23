#' This function validates the parameters for RRW.
#'
#' This is an internal function that validates the parameters for RRW.  It is used with assessRRWfit(), which is used with an opimization program.  Validate parameters ensures that the parameters input are valid. There are no default values because all values must be assessed for validity. There are no defaults, so all the parameters can be checked for validity.
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param b A vector of number(s) specifying the boundary distance from a 0 startpoint. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bCols" argument.
#' @param bcs A vector of positive values that specify the likelihood that the boundary will be reduced because there it is taking too many samples before reaching a threshold. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of b. The column names must be specified in the "bcsCols" argument. Essentially, this models a decision-maker reducing the amount of information that they need before a response because it is taking to long.  The boundary reduces when there is little change in the likelihood of a response after N number of steps.  The routine makes this determination as a function of the size of the unaltered boundary (boundaryChangeSensitivity * boundary). Small values (e.g., 0.1), result in an impatient person (relatively early chnage in the boundary values).  Large values (e.g., 0.9), result in a patient person (a boundary that is relatively resistent to change).  A value of 0 is an exception: it will result in no boundary reduction under any circumstances. Impatience increases the likelihood of an error, so it increases the average time of error responses.
#' @param Ter A vector of positive values that adjusts the time of encoding/response (non-decision time). When Ter = 0, a single, best fitting Ter will be calculated. Use this parameter when you have several conditions and you hypothesize that some conditions should have a larger Ter than others.  For example, this may model longer encoding process when different numbers of items are on the screen for different conditions. The hypothetical fastest condition should have Ter = 0, which simply indicates that the RRW program will find the best fit Ter. The slower conditions will have larger Ter values and will be fit ralative to the fastest condition. The change in samples resulting from the Ter effect is a function of the size of the boundary (Ter * boundary). This is still experimental.
#' @param s A vector of signed number(s) between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.   If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of s. The column names must be specified in the "sCols" argument.
#' @param nSD A vector of positive number(s) representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step.  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD. The column names must be specified in the "nSDCols" argument.
#' @param db A vector of signed number(s) representing the beta value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of db. The column names must be specified in the "dbCols" argument.
#' @param da A vector of positive number(s) representing the asymptote value in the Information Accrual Bias (IAB).  If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of da. The column names must be specified in the "daCols" argument.
#' @param vc A vector of signed number(s) that indicates the change of overlap that one predicts. If the number of values is greater than 1, then each value must have a corresponding effect coded column in the dataset. These columns should be effect coded, because the values contained in this column will be multiplied by the value of vc. The column names must be specified in the "vcCols" argument.
#' @param bCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of boundary (b).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of b.
#' @param bcsCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the boundaryChangeSensitivity (bcs).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of bcs. All effect/dummy codes must be positive.
#' @param TerCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the non-decision time (Ter).  These columns should be effect/dummy coded, because the values contained in this column will be multiplied by the value of Ter.  All effect/dummy codes must be positive.
#' @param sCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of startValue (s).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of s.
#' @param nSDCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of noiseSD (nSD).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of nSD.
#' @param dbCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decayBeta (db).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of db.
#' @param daCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of decay asymptote (da).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of da.
#' @param vcCols A vector of strings that identifies the name(s) of the column in data that identifies the conditions that will affect the change of value change (vc).  These columns should be effect coded, because the values contained in this column will be multiplied by the value of overlap.
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @return TRUE if the parameters are valid; FALSE if the parameters are invalid
#' @keywords RRW random walk parameters validate
#' @export
#' @examples validateParameters (data, b=14, s=0.1, nSD = 0, db = .2, da = .2, vc = 0, sCols=NULL, nSDCols=NULL, dbCols=NULL, daCols=NULL, vcCols="valueChangeEffect")

validateParameters <- function(data, b, bcs, Ter, s, nSD, db, da, vc, bCols, bcsCols, TerCols, sCols, nSDCols, dbCols, daCols, vcCols, dataOverlapCol) {

  out <- TRUE

  dataRunCols <- c(dataOverlapCol, bCols, sCols, nSDCols, dbCols, daCols, vcCols, bcsCols, TerCols)

  #extract the unique information
  data.tmp <- unique(data[,dataRunCols, drop=F])
  data.n <- nrow(data.tmp)

  ovIn <- NULL
  bIn <- NULL
  bcsIn <- NULL
  bcsIn <- NULL
  TerIn <- NULL
  sIn <- NULL
  nSDIn <- NULL
  daIn <- NULL
  dbIn <- NULL

  for(i in 1:data.n) {
    ovIn[i] <- data.tmp[i,dataOverlapCol] + getRowParameterValue(data.tmp[i,], vcCols, vc)
    bIn[i] <- getRowParameterValue(data.tmp[i,], bCols, b)
    bcsIn[i] <- getRowParameterValue(data.tmp[i,], bcsCols, bcs)
    TerIn <- getRowParameterValue(data.tmp[i,], TerCols, Ter)
    sIn[i] <- getRowParameterValue(data.tmp[i,], sCols, s)
    nSDIn[i] <- getRowParameterValue(data.tmp[i,], nSDCols, nSD)
    daIn[i] <- getRowParameterValue(data.tmp[i,], daCols, da)
    dbIn[i] <- getRowParameterValue(data.tmp[i,], dbCols, db)
  }

  if(min(bIn) <= 1) out <- FALSE
  if(min(bcsIn) < 0) out <- FALSE
  if(min(TerIn) < 0) out <- FALSE
  if(max(ovIn) > 2) out <- FALSE
  if(max(abs(sIn)) >= 1) out <- FALSE
  if(min(nSDIn) < 0) out <- FALSE

  return(out)

}
