#' Function that standardizes the data and fit of the simulation
#'
#' This function standardizes the data and fit of the simulation
#' @param data a vector of numbers that are the behavioral data.
#' @param fit a vector of numbers that were fit to the data by a model.
#' @return a dataframe with the standardized data and fit
#' @keywords standardize scale
#' @export
#' @examples standardizeDataAndFit (data, fit)


standardizeDataAndFit <- function(data, fit) {

  data.all <- c(data, fit)
  data.all.m <- mean(data.all, na.rm = T)
  data.all.sd <- sd(data.all, na.rm = T)
  data.z <- (data - data.all.m)/data.all.sd
  sim.z <- (fit - data.all.m)/data.all.sd
	
	df.out <- data.frame(data = data.z, fit = sim.z)
	
	return(df.out)

}