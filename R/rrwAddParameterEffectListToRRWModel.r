
#' This  function is used to create or modify an rrwModelList
#'
#' Function that creates or modifys an rrwModelList by adding effects to the list one at a time.
#' @param rrwModelList A list that specifies the rrw model.  If you are building a new list, then this should be NULL. If you are adding a new effect to an old list, then this should be the rrwModelList that you are adding an effect to. DEFAULT = NULL
#' @param parameterEffectList A list that specifies the parameter effect (or a list of effects, e.g., c(effect1, effect2, effect3)).  DEFAULT = NULL
#''
#' @return the updated rrwModelList
#' @keywords rrw rrwModelList update add effect parameter effect
#' @export
#' @examples rrwAddParameterEffectListToRRWModel (modelList, parameterEffectList)

rrwAddParameterEffectListToRRWModel <- function(rrwModelList = NULL, parameterEffectList = NULL) {

  parameter <- names(parameterEffectList)
  for(i in 1:length(parameter)) {
    #if the rrw modelList is NULL, create it.
    if(is.null(rrwModelList)) {
      rrwModelList <- parameterEffectList[i]
    } else {
      #if the rrw modelList is not NULL, but it does not yet contain this parameter.
      if(is.null(rrwModelList[[parameter[i]]])) {
        rrwModelList[[parameter[i]]] <- parameterEffectList[[i]]
      } else {
        rrwModelList[[parameter[i]]] <- c(rrwModelList[[parameter[i]]], parameterEffectList[[i]])
      }
    }
  }

    return(rrwModelList)
}
