#' This function modifies an rrwModelList so that the parameters values are fixed to the output of the parameter values in a runStats list (which is the output of the rrw fit functions). It is used to test a fixed parameter model on data using the model parameters extracted from another dataset.
#'
#' Function that modifies the parameter values of an rrwModelList to those of a runStats list.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param runStats A list that contains the output of a single run of the kFold cross validation of the RRW.  It is output by rrw fit functions and used in the rrwKfoldStatsToDataframe, etc.
#''
#' @return a rrwModelList that has the parameters values are fixed to the output of the parameter values in a runStats list.
#' @keywords rrw RRW modelList runStats
#' @export
#' @examples rrwConvertFreeRRWModelToFixedRRWModel (rrwModelList,runStats)

rrwConvertFreeRRWModelToFixedRRWModel <- function(rrwModelList, runStatsList) {
  fixedParameterList <- runStatsList$parameters
  freeModelParameters <- names(rrwModelList)

  for(i in freeModelParameters) {
    nParamsFree <- length(rrwModelList[i])
    nParamsFix <- length(fixedParameterList[i])
    for (j in 1:nParamsFree) {
      #get the column name of each instance of this parameter
      freeColName <- rrwModelList[[i]][[j]][["columnName"]]
      #find the same column name in the fixed model
      if(is.null(freeColName)) {
        freeColName <- i
      }
      for(k in 1:nParamsFix) {
        if (fixedParameterList[[i]][["columns"]][[k]] == freeColName) {
          rrwModelList[[i]][[j]][["parameterHigh"]] <- fixedParameterList[[i]][["values"]][[k]]
          rrwModelList[[i]][[j]][["parameterLow"]] <- fixedParameterList[[i]][["values"]][[k]]
        }
      }
    }
  }
  return(rrwModelList)
}
