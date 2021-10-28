#' This function uses the rrwModelList to collapse the raw data into usable RRW data
#'
#' Function that uses the rrwModelList to collapse the raw data into usable RRW data
#' @param data This is a dataframe containing the raw trial-by-trial data that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param grpCols A vector of strings containing the grouping columns. If NULL this argument will be ignored.  DEFAULT = NULL
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param correctCol a string that specifies the name of the column that specifies if the participant chose the item with the greatest value distribution (correct) or if they did not (incorrect).
#' @param correctVals a vector of two values that specifies the "correct" value (index 1) and the "incorrect" value (index 2). e.g, c("yes", "no"). DEFAULT = c(TRUE, FALSE)
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param splitByCorrect A boolean that specifies whether to use the "correct/incorrect" column as a grouping variable for the RT data. DEFAULT = TRUE.
#''
#' @return The input data for the RRW, which is, essentially, the data summarized by the critical condition, overlap, etc.
#' @keywords RRW format data input
#' @export
#' @importFrom dplyr %>%
#' @examples rrwRawDataToRRWFormatData (data, modelList, "rt", "correct", c(TRUE, FALSE), "overlapRound")

rrwRawDataToRRWFormatData <- function(data, grpCols = NULL, dataRtCol = "rt", correctCol = "correct",correctVals = c(TRUE, FALSE), dataOverlapCol = "overlapRound", splitByCorrect = TRUE) {

  #ensure correct/incorrect is coded as 1/0
  data$correct01 <- ifelse (data[[correctCol]]==correctVals[1], 1, 0)

  ddplyGroupColumns <- c(dataOverlapCol, grpCols)

  #here are the columns we will keep in the end.  We add it here, because later
  #it will be more difficult. The columns in "quotes" are columns we create inside this function
  keepCols <- c(ddplyGroupColumns, "pHit", "rt", "n")

  #here we collapse the data for the correct and incorrect and RT data.
  #correct == T
  df.rwPhit.T<-data.frame(data %>% dplyr::group_by_at(ddplyGroupColumns) %>% dplyr::summarise(pHit = mean(correct01, na.rm=T)))
  #correct == F
  df.rwPhit.F<-data.frame(data %>% dplyr::group_by_at(ddplyGroupColumns) %>% dplyr::summarise(pHit = 1 - mean(correct01, na.rm=T)))
  #prepare RT data
  if(splitByCorrect) {
    ddplyGroupColumns <- c(ddplyGroupColumns, "correct01")
    keepCols <- c(keepCols, "correct")
  }
  df.rwRT<-data.frame(data %>% dplyr::group_by_at(ddplyGroupColumns) %>% dplyr::summarise(rt = mean(eval(parse(text = dataRtCol)), na.rm=T), n = length(correct01)))

  if(splitByCorrect) {
    #we now code the correct and incorrect using TRUE and FALSE
    df.rwPhit.T$correct <- TRUE
    df.rwPhit.F$correct <- FALSE
    df.rwRT$correct <- ifelse (df.rwRT$correct01 == 1, TRUE, FALSE)

    #merge correct and incorrect data for phit dataset
    df.pHit <- merge(df.rwPhit.T,df.rwPhit.F, all=T)
  } else {
    df.pHit <- df.rwPhit.T
  }

  #merge pHit and RT data
  df.raw <- merge(df.pHit,df.rwRT, all=T)

  #keep only needed columns
  df.raw <- df.raw[, keepCols]

  #return the new, collasped RRW data
  return(df.raw)
}
