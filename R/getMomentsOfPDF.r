#' This function fits returns the moments of a Probability Distribution Function.
#'
#' This function samples a Probability Distribution Function, and returns the moments of that PDF.
#' @param PDFps A vector containing the proportions of the PDF.
#' @param values A vector containing the values (x-axis) of the CDF.
#' @param samples A number containing the number of samples to be used in the simulation. DEFAULT = 50000.
#' @return A dataframe that contains the moments of the simulated values from the PDF ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile)
#' @keywords PDF probability distribution function moments values
#' @export
#' @examples getMomentsOfPDF (PDFvector, valueVector)

getMomentsOfPDF <- function(PDFps, values, samples = 50000) {
  tmpDat <- as.numeric(sample(values, samples, prob=PDFps, replace=T))
  df.out <- data.frame(mean = mean(tmpDat), sd = sd(tmpDat), Q25 = as.numeric(quantile(tmpDat, .25)), Q50 = median(tmpDat), Q75 = as.numeric(quantile(tmpDat, .75)))
  return(df.out)

}
