#' This function extracts the parameter bounds, parameter column names, and number of free parameters from the rrwModelList.
#'
#' Function that extracts the parameter bounds, parameter column names, and number of free parameters from the rrwModelList.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param startPars.n An integer that specifies the number of free parameters that are not specified in the modelList. This will be added to the total number of free parameters in the modelList.  DEFAULT = 1 (#Start at 1 because of Ter)
#''
#' @return a list containing the parameter column names (pCols), the number of free parameters (pars.n), and the parameter upperBound, lowerBound, and interval (pU, pL, pI)
#' @keywords rrw RRW parameter specs extract
#' @export
#' @examples rrwGetModelParameterSpecs (modelList)

rrwGetModelParameterSpecs <- function (rrwModelList, startPars.n = 1) {
  #Specify all the possible parameters
  #parameters <- c("s", "b", "db", "da", "nSD", "vc")

  #Get the parameters from the modelList
  parameters <- names(rrwModelList)

  #Initiate all the lists that will be output
  pU <- list()
  pL <- list()
  pI <- list()
  pCols <- list()
  #initiate the number of free parameters
  pars.n <- 1 #Start at 1 because of Ter

  #for each paramter
  for(i in parameters) {
    #is there a specified parameter in the ModelList?
    if(!is.null(rrwModelList[[i]])) {
      xU <- NULL
      xL <- NULL
      xI <- NULL
      xCol <- NULL
      #if so, for each sublist do the following
      for(j in 1:length(rrwModelList[[i]])) {
        #Extract the high, low, and interval bounds for the parameters
        tmpH <- rrwModelList[[i]][[j]]$parameterHigh
        #and give them the appropriate ParameterName
        names(tmpH) <- rrwModelList[[i]][[j]]$ParameterName
        xU <- c(xU, tmpH)
        tmpL <- rrwModelList[[i]][[j]]$parameterLow
        names(tmpL) <- rrwModelList[[i]][[j]]$ParameterName
        xL <- c(xL, tmpL)
        tmpI <- rrwModelList[[i]][[j]]$parameterInt
        names(tmpI) <- rrwModelList[[i]][[j]]$ParameterName
        xI <- c(xI, tmpI)

        #if the high > low, then it is a free parameter, so add 1 to pars.n
        pars.n <- ifelse (tmpH > tmpL, pars.n + 1, pars.n)

        #if there is a column in the rrwData, add it to a column variable here
        if(!is.null(rrwModelList[[i]][[j]]$columnName)) {
          tmp <- rrwModelList[[i]][[j]]$columnName
          names(tmp) <- rrwModelList[[i]][[j]]$ParameterName
          xCol <- c(xCol,tmp)
        }
      }
      if(!is.null(xCol)) {
        pColN <- paste(i,"Cols", sep="")
        pCols[[pColN]] <- xCol
      }
      pU[[i]] <- xU
      pL[[i]] <- xL
      pI[[i]] <- xI
    }
  }
  names(pars.n) <- "pars.n"

  #create an outlist
  outList <- list(pCols = pCols, pars.n = pars.n, pU = pU, pL = pL, pI = pI)

  #return it
  return(outList)
}
