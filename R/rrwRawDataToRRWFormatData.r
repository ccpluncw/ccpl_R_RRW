#' This function uses the rrwModelList to collapse the raw data into usable RRW data
#'
#' Function that uses the rrwModelList to collapse the raw data into usable RRW data
#' @param data This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param correctCol a string that specifies the name of the column that specifies if the participant chose the item with the greatest value distribution (correct) or if they did not (incorrect).
#' @param correctVals a vector of two values that specifies the "correct" value (index 1) and the "incorrect" value (index 2). e.g, c("yes", "no"). DEFAULT = c(TRUE, FALSE)
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#''
#' @return The input data for the RRW, which is, essentially, the data summarized by the critical condition, overlap, etc.
#' @keywords RRW format data input
#' @export
#' @examples rrwRawDataToRRWFormatData (data, modelList, "rt", "correct", c(TRUE, FALSE), "overlapRound")

rrwRawDataToRRWFormatData <- function(data, rrwModelList, dataRtCol = "rt", correctCol = "correct",correctVals = c(TRUE, FALSE), dataOverlapCol = "overlapRound") {

  #ensure correct/incorrect is coded as 1/0
  data$correct01 <- ifelse (data[[correctCol]]==correctVals[1], 1, 0)

  #identify the group_by collumns for the ddply
  #start by putting in the overlap column
  ddplyGroupColumns <- dataOverlapCol
  #add the grouping columns specified in the rrwModelList
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
  ddplyGroupColumns <- c(ddplyGroupColumns, grpCols)

  #here are the columns we will keep in the end.  We add it here, because later
  #it will be more difficult. The columns in "quotes" are columns we create inside this function
  keepCols <- c(ddplyGroupColumns, "correct", "pHit", "rt", "n")

  #here we collapse the data for the correct and incorrect and RT data.
  #correct == T
  df.rwPhit.T<-data.frame(data %>% group_by_at(ddplyGroupColumns) %>% summarise(pHit = mean(correct01)))
  #correct == F
  df.rwPhit.F<-data.frame(data %>% group_by_at(ddplyGroupColumns) %>% summarise(pHit = 1 - mean(correct01)))
  #prepare RT data
  ddplyGroupColumns <- c(ddplyGroupColumns, "correct01")
  df.rwRT<-data.frame(data %>% group_by_at(ddplyGroupColumns) %>% summarise(rt = mean(!!sym(dataRtCol)), n = length(correct01)))

  #we now code the correct and incorrect using TRUE and FALSE
  df.rwPhit.T$correct <- TRUE
  df.rwPhit.F$correct <- FALSE
  df.rwRT$correct <- ifelse (df.rwRT$correct01 == 1, TRUE, FALSE)

  #merge correct and incorrect data for phit dataset
  df.pHit <- merge(df.rwPhit.T,df.rwPhit.F, all=T)

  #merge pHit and RT data
  df.raw <- merge(df.pHit,df.rwRT, all=T)

  #keep only needed columns
  df.raw <- df.raw[, keepCols]

  #return the new, collasped RRW data
  return(df.raw)
}
