#' This  function uses the rrwModelList to add the effect and dummy variables to the RRW data (one column at a time)
#'
#' Function that uuses the rrwModelList to add the effect and dummy variables to the RRW data (one column at a time)
#' @param data This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#''
#' @return the data is returned with the added dummy/effect coded column
#' @keywords rrw RRW data effect dummy code column
#' @export
#' @examples rrwCodeColumn (data, modelList)

rrwCodeColumn <- function(data,rrwModelList) {

  #first identify the number of conditions in the new coded variable (the number of rows in df.code)
  numConds <- length(rrwModelList$df.code[,1])

  #make sure the last condition is "default"  It is a necessary condition
  if(rrwModelList$df.code[numConds,1] != 'default') {
    stop("last df.code entry must be 'default' with a value")
  }

  #if there is only one condition (i.e., default)
  if(numConds == 1) {
    #code the column as a constant equal to the default
    data[[rrwModelList$columnName]] <- rrwModelList$df.code[numConds,2]
  } else {

    #code the first condition, make all the rest equal to the default
    data[[rrwModelList$columnName]] <- with(data,ifelse(eval(parse(text = as.character(rrwModelList$df.code[1,1]))), rrwModelList$df.code[1,2], rrwModelList$df.code[numConds,2]))
    #if there are more than 2 conditions
    if(numConds > 2) {
      #now, change each other condition from the default to the appropriate code (one at a time)
      for(i in 2:(numConds -1)) {
        data[[rrwModelList$columnName]] <- with(data,ifelse(eval(parse(text = as.character(rrwModelList$df.code[i,1]))), rrwModelList$df.code[i,2], data[[rrwModelList$columnName]]))
      }
    }
  }
  #return the coded data with the added column
  return(data)
}
