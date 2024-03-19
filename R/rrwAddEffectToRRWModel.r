
#' This  function is used to create or modify an rrwModelList
#'
#' Function that creates or modifys an rrwModelList by adding effects to the list one at a time.
#' @param rrwModelList A list that specifies the rrw model.  If you are building a new list, then this should be NULL. If you are adding a new effect to an old list, then this should be the rrwModelList that you are adding an effect to. DEFAULT = NULL
#' @param parameter A string that specifies the parameter that the effect being added will influence.  It can take on one of the following values: "s","b","nSD", "db". "da". "vc" DEFAULT = NULL
#' @param columnName A string that specifies the name of the column that will be created and added to the data.  The column will contain the effect/dummy code of this effect. This columnName can remain NULL if this effect is constant across all conditions (therefore the dummycode would equal 1 for all conditions and therefore would be redundant). DEFAULT = NULL
#' @param ParameterName A string that specifies the name of the free parameter for this effect.  It is not the "parameter" argument above.  It is, rather, a variable name that holds the free parameter values and is output in the output file. You will identify this influence of this effect by the ParameterName. DEFAULT = NULL.
#' @param parameterBounds A vector that specifies the upperBound, lowerBound, and minimum interval, c(upper, lower, interval), that the free parameter will vary when optimized. If the upperBound==lowerBound, then this will be a fixed parameter, that is fixed at the upperBound value. DEFAULT = c(0,0,.001)
#' @param df.code A dataframe that specifies the conditional statements for coding the dummy or effect variable.  The first column is a string stating the conditional (e.g., "variableName == 'level'"), the second column is a numeric stating the value (1). Each row, is a different conditional/value. The last row, first column, must be "default" and the second column must be the default value.  This df.code can remain NULL if this effect is constant across all conditions (therefore there are no conditionals). DEFAULT = NULL.
#' @param GroupByVariables A vector of strings that specify the variable names of the grouping variable(s) that will be used by the ddply to create the summary dataset.  DEFAULT = NULL.
#''
#' @return the updated rrwModelList
#' @keywords rrw rrwModelList update add effect
#' @export
#' @examples rrwAddEffectToRRWModel (modelList, parameter = "b", columnName = NULL, ParameterName = overallB, parameterBounds = c(100,10,0.001), df.code = NULL, GroupByVariables = NULL))

rrwAddEffectToRRWModel <- function(rrwModelList = NULL, parameter = NULL, columnName = NULL, ParameterName = NULL, parameterBounds = c(0,0,0.001), df.code = NULL, GroupByVariables = NULL) {

  validParameters <- c("s", "b", "bcs", "Ter", "db", "da", "nSD", "vc", "ec", "rdN", "rdMean", "rdSD")

  if(!(parameter %in% validParameters)) {
    stop(paste("parameter option must take on one of the following values: ", paste(validParameters, sep="", collapse=" ")))
  }

  if( (!is.null(df.code) & is.null(columnName)) ) {
    stop("columnName option must have a value (string) because you are specifying the values of that column in the df.code option. columnName option is not necessary if this parameter is a fixed value.  In that case, then df.code is not necessary either.")
  }

  if( parameterBounds[1] < parameterBounds[2] ) {
    stop("parameterBounds[1] (the high boundary) must be equal to or greater than parameterBounds[2] (the low boundary). If parameterBounds[1] (the high boundary) is equal to parameterBounds[2] (the low boundary), then this parameter will be fixed at the value of parameterBounds[1]")
  }

  if( parameterBounds[3] <= 0 ) {
    stop("parameterBounds[3] (the interval) must be greater than 0 even if this parameter is fixed.  If it is fixed, then set parameterBounds[1] (the high boundary) equal to parameterBounds[2] (the low boundary).")
  }

  if( (!is.null(df.code) & is.null(GroupByVariables)) ) {
    print("Note: GroupByVariables option has no value. If df.code is conditioning its values on the levels of a particular variable or variables, they should be listed in the GroupByVariables option.  If df.code is simply a constant (unconditioned) value, then the GroupByVariables option is not needed.")
  }

  if(!is.null(df.code)) {
    numConds <- length(df.code[,1])
    #make sure the last condition is "default"  It is a necessary condition
    if(df.code[numConds,1] != 'default') {
      stop("last df.code entry must be 'default' with a value")
    }
  }

  tmpList <- list (columnName = columnName, df.code = df.code, GroupByVariables = GroupByVariables, ParameterName = ParameterName, parameterHigh = parameterBounds[1], parameterLow = parameterBounds[2], parameterInt = parameterBounds[3])

  if(is.null(rrwModelList[[parameter]])) {
    rrwModelList[[parameter]] <- list(tmpList)
  } else {
    rrwModelList[[parameter]] <- c(rrwModelList[[parameter]], list(tmpList))
  }

  return(rrwModelList)
}
