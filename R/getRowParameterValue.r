#' This function outputs the parameter value for a specific row of the data.
#'
#' This function outputs the parameter value for a specific row of the data. Because models are effect coded, each parameter value for the data row is often the sum of the effects.  This function calculates and outputs that sum.
#' @param data.row This is a single row of the  dataframe that  contains (at least) the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param parCol This is a vector of strings that identify the names of the columns in data that contain the effect codes for the specific parameter entered here. Default  = NULL (e.g., no effect codes for this parameter)
#' @param parVec A vector of number(s) specifying the values for this specific paramter.
#' @return The value(s) of the paramter for the condition specified in data.row.
#' @keywords RRW random walk parameters effect codes
#' @export
#' @examples getRowParameterValue (data[1,], "boundaryCols", b)

getRowParameterValue <- function(data.row, parCol=NULL, parVec) {

  if(!is.null(parCol)) {
    if(length(parVec) != length(parCol)) {
      print("parameters: ")
      print(parVec)
      print(" columns: ")
      print(parCol)
      stop("number of parameters must equal the number of columns in the data devoted to those parameters")
    }

    #add all the effects for this parameter and return
    tmp.par <- 0
    for (j in 1:length(parCol)) {
      tmp.par <- tmp.par + data.row[[parCol[j]]] * parVec[j]
    }
  } else {
    if(length(parVec) != 1) {
      print("parameters: ")
      print(parVec)
      stop("number of parameters must equal 1 if there are no columns in the data devoted to those parameters")
    }
    #if no effects, then it is a single parameter value that simply varies
    tmp.par <- parVec
  }

  return(tmp.par)

}
