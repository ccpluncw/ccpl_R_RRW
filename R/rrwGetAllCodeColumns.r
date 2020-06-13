#' This  function uses the rrwModelList to add all the effect and dummy variables to the RRW data
#'
#' Function that uses the rrwModelList to add all the effect and dummy variables to the RRW data
#' @param data This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#''
#' @return the data is returned with all the added dummy/effect coded columns
#' @keywords rrw RRW data effect dummy code column
#' @export
#' @examples rrwGetAllCodeColumns (data, modelList)

rrwGetAllCodeColumns <- function(data,rrwModelList) {

  #add all the coded columns to the dataset

  #For each parameter in the ModelList,
    for(i in 1:length(rrwModelList)) {
      #for each sublist in a parameter
      for(j in 1:length(rrwModelList[[i]])) {
        #for each column
        if(!is.null(rrwModelList[[i]][[j]]$columnName)) {
          #for each parameter: tmpStartModel [[i]]
            #for each column in the parameter: tmpStartModel [[-]] [[j]]
          data <- rrwCodeColumn(data, rrwModelList[[i]][[j]])
        }
      }
    }
  return(data)
}
