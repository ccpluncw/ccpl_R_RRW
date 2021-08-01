#' A function to get the grouping columns from the rrwModelList.
#'
#' This function gets the grouping columns from the rrwModelList.
#'
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @keywords morals RT pHit fit plots
#' @return a vector of strings that are the names of the grouping column headers.
#' @export
#' @examples getGroupColsFromModelList (simpleModelList)

getGroupColsFromModelList <- function(rrwModelList) {
  grpCols <- NULL
  #try to do it with nested lists
  tryCatch ({
    grpCols <- unique(as.vector(unlist(lapply(unlist(rrwModelList, recursive=FALSE), `[`, "GroupByVariables"))))
    }, error = function(e) {}
  )
  #if that fails, try unnested lists
  if(is.null(grpCols)) {
    tryCatch ({
      grpCols <- unique(as.vector(unlist(lapply(rrwModelList, `[`, "GroupByVariables"))))
      }, error = function(e) {}
    )
  }
  return(grpCols)
}
