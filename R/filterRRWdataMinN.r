#' This function uses the rrwModelList to collapse the raw data into usable RRW data
#'
#' Function that uses the rrwModelList to collapse the raw data into usable RRW data
#' @param data This is a dataframe that is used as the input data for the RRW, which is, essentially, the data summarized by the critical condition, overlap, etc. It is output by rrwRawDataToRRWFormatData().
#' @param grpCols A vector of strings containing the grouping columns. If NULL this argument will be ignored.  DEFAULT = NULL
#' @param dataRtCol A string that identifies the name of the column in data that contains the average RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials for the specific overlap/correct/condition combination. Default is "pHit"
#' @param nCol A string that identifies the name of the column in data that contains the number of trials for the specific overlap/correct/condition combination. Default is "n"
#' @param dataOverlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlapRound"
#' @param minN An integer that specifies the minimum number of trials required in the calculation of summary statistic in the dataRtCol and dataPhitCol for the analysis to run on that row of data. DEFAULT = 20.
#''
#' @return The input data for the RRW, filtered by minN
#' @keywords RRW format data input
#' @export
#' @importFrom dplyr %>%
#' @examples filterRRWdataMinN (data, c("condition"), "rt", "pHit", "n", "overlapRound", 40)

filterRRWdataMinN <- function(data, grpCols = NULL, dataRtCol = "rt", dataPhitCol = "pHit", nCol = "n", dataOverlapCol = "overlapRound", minN = 20) {

  data[[dataRtCol]] <- ifelse(data[[nCol]] < minN, NA, data[[dataRtCol]])

  ddplyGroupColumns <- c(dataOverlapCol, grpCols)
  df.tmp.pHit.sum <- data %>% group_by_at(ddplyGroupColumns) %>% summarize(N = sum(eval(parse(text = nCol)), na.rm=T))
  df.tmp.pHit.sum$keepPhit <- ifelse(df.tmp.pHit.sum$N > minN, TRUE, FALSE)
  data <- merge(data,df.tmp.pHit.sum, all=T )
  data[[dataPhitCol]] <- ifelse(data$keepPhit, data[[dataPhitCol]], NA)
  data$keepPhit <- NULL
  data$N <- NULL

  return(data)
}
