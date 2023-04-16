#' This function converts the runStats list (which is the output of the rrw fit functions) into a single row of a dataframe.
#'
#' Function that converts the runStats list into a single row of a dataframe.
#' @param runStats A list that contains the output of a single run of the kFold cross validation of the RRW.  It is output by rrw fit functions and used in the rrwKfoldStatsToDataframe, etc.
#''
#' @return a single row of a dataframe that contains the parameter values and fit statistics in each column.
#' @keywords rrw RRW kFold runStats
#' @export
#' @examples rrwKfoldStatsToDataframe (runStats)


rrwRunStatsToDataframe <- function (runStats) {

  df.boundary.1 <- runStats$parameters$b
  df.boundary.2 <- data.frame(t(df.boundary.1$values))
  colnames(df.boundary.2) <- df.boundary.1$columns

  df.StartValue.1 <- runStats$parameters$s
  df.StartValue.2 <- data.frame(t(df.StartValue.1$values))
  colnames(df.StartValue.2) <- df.StartValue.1$columns

  df.NoiseSD.1 <- runStats$parameters$nSD
  df.NoiseSD.2 <- data.frame(t(df.NoiseSD.1$values))
  colnames(df.NoiseSD.2) <- df.NoiseSD.1$columns

  df.DecayBeta.1 <- runStats$parameters$db
  df.DecayBeta.2 <- data.frame(t(df.DecayBeta.1$values))
  colnames(df.DecayBeta.2) <- df.DecayBeta.1$columns

  df.DecayAsymptote.1 <- runStats$parameters$da
  df.DecayAsymptote.2 <- data.frame(t(df.DecayAsymptote.1$values))
  colnames(df.DecayAsymptote.2) <- df.DecayAsymptote.1$columns

  df.ValueChange.1 <- runStats$parameters$vc
  df.ValueChange.2 <- data.frame(t(df.ValueChange.1$values))
  colnames(df.ValueChange.2) <- df.ValueChange.1$columns

  df.ter.1 <- data.frame(ter = runStats$parameters$ter$ter)
  df.AIC.1 <- data.frame(AIC = runStats$fitStats$AIC)
  df.BIC.1 <- data.frame(BIC = runStats$fitStats$BIC)
  df.r2.1 <- data.frame(r2 = runStats$fitStats$r2)

  df.tmp.stats <- data.frame(df.boundary.2, df.StartValue.2, df.NoiseSD.2, df.DecayBeta.2, df.DecayAsymptote.2, df.ValueChange.2, df.ter.1, df.AIC.1, df.BIC.1, df.r2.1)

  return(df.tmp.stats)
}
