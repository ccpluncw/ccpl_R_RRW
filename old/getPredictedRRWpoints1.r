#' This function simulates the RRW for a single set of parameters across all values of overlap.
#'
#' This function simulates the RRW for a single set of parameters across all values of overlap and then attaches this simulation to the empirical RT and error data.  This simulation's fit will be assessed by the optimization program (using assessRRWfit()) and then the parameters will be adjested.
#'
#'
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contain a column specifying a condition that will influence either the startpoint or the value of the trials.
#' @param b A number specifying the boundary distance from a 0 startpoint. This value is specific to the RRW simulation and has no default value
#' @param RWkeepColumns A vector of strings that specify the columns of the RRW simulated data that should be kept when the simulated dataframe is merged with the empirical dataframe.
#' @param mergeByDataColumns A vector of strings that specify the columns of the empirical data that should be used in the "by.x" option of the merge() function when the simulated dataframe is merged with the empirical dataframe.
#' @param mergeByRWColumns A vector of strings that specify the columns of the RRW simulated data that should be used in the "by.y" option of the merge() function when the simulated dataframe is merged with the empirical dataframe.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param RwSamplesCol A string that identifies the name of the column in the RRW simulations that contains the summary of the samples that you want to use as a simulation for RT.  The possible columns are: "Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile).  The default is "Q50" because the median is more robust than the mean.
#' @param startValue A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  startValue = 0 is the default and represents an unbiased start point.
#' @param noiseSD A positive number representing the SD of the noise distribution.  The noise distribution is N(0,nSD) and is added to the value of every step. noiseSD = 0 is the default and represents no noise being added to each step.
#' @param decayAsymptote A positive number representing the asymptote value in the Information Accrual Bias (IAB). decayAsymptote = 0.2 is the default.
#' @param decayBeta A signed number representing the beta value in the Information Accrual Bias (IAB). decayBeta = 0 is the default and represents no IAB.
#' @param startBiasEffect A signed number between -1 and 1 that indicates the position of the start point as a proportion of the boundary value.  The startBiasEffect is a constant that is added to s (the start point) as a function of the values contained in the startBiasEffectCol. The default is NULL because there is no startBiasEffect unless one includes a startBiasEffectCol.
#' @param startBiasEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the position of the start point.  This column should be effect coded, because the values contained in this column will be multiplied by the value of startBiasEffect (and then added to s) to get the start point of the RRW.
#' @param valueChangeEffect A signed number between that indicates the change of overlap that one predicts as a function of the values contained in the valueChangeEffectCol. The default is NULL because there is no valueChangeEffect unless one includes a valueChangeEffectCol.
#' @param valueChangeEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the change of overlap (value change).  This column should be effect coded, because the values contained in this column will be multiplied by the value of valueChangeEffect (and then added to overlap) to get the overlap input in the RRW. This change is assumed to be constant across the entire overlap sequence (from 0-1).
#' @param loops A number specifying the number of loops that will be run in the RRW simulation when it calculates the summary statistics for each number of samples for each boundary. Higher numbers produce more precise estimates, but also increase the time needed to converge on a solution.  Default is 200.
#' @param progress TRUE or FALSE that specifies whether to present a progress bar.  Default is FALSE.
#'
#' @return A dataframe that contains the "data" plus the fitted values from the model ("Q25" (the 25th quartile); "Q50" (the median); "mean" (the mean); "Q75" (the 75th quartile); pCross (the fitted pHit from the model - probability of crossing each boundary); rtFit (the fitted rt values from the model))
#' @keywords RRW random walk simulation get predicted
#' @export
#' @examples getPredictedRRWpoints (data, b=14, RWkeepColumns = c("overlap", "Q50", "pCross", "correct"), mergeByDataColumns = c("overlap", "correct"), dataRtCol = "rt", RwSamplesCol = "Q50", s=0.1, loops = 400)

getPredictedRRWpoints1 <- function (data, b, RWkeepColumns, mergeByDataColumns, mergeByRWColumns, dataRtCol, RwSamplesCol, startValue = 0,  noiseSD = 0, decayAsymptote = 0.2, decayBeta = 0.0, startBiasEffect = NULL, startBiasEffectCol = NULL, valueChangeEffect = NULL, valueChangeEffectCol = NULL, loops = 200, progress = F) {

    overlapSeq <- unique(data$overlap)


    if(!is.null(startBiasEffectCol) & !is.null(valueChangeEffectCol)) {
      #if there is a start bias effect parameter, get the start bias effect conditions.  This should be either dummy coded or effect coded.
      sbVcConds <- unique(data[c(startBiasEffectCol, valueChangeEffectCol)])
      sbVcConds.n <- nrow(sbVcConds)

      #now, for each start bias effect X value change condition, run the simple random walk with the appropriate parameters and rbind them together
      df.out <- NULL
      for(i in 1:sbVcConds.n) {
        sC <- sbVcConds[i,startBiasEffectCol]
        vC <- sbVcConds[i,valueChangeEffectCol]

        sIn <- startValue + (sC * startBiasEffect)
        overlapSeqIn <- overlapSeq + (vC * valueChangeEffect)

        df.tmp <- getMomentsOfRRWoverlap(overlapSeqIn, b, startValue = sIn, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops = loops, progress=progress)
        df.tmp[[startBiasEffectCol]] <- sC
        df.tmp[[valueChangeEffectCol]] <- vC

        #replace distorted overlaps (overlapSeq + vCp*valueChangeEffect) with the original overlaps, so it looks like real data
        df.ov <- data.frame (overlap = overlapSeqIn, origOverlap = overlapSeq)
        df.tmp <- merge (df.tmp, df.ov, by="overlap")
        colnames(df.tmp)[colnames(df.tmp)=="overlap"] <- "distOv"
        colnames(df.tmp)[colnames(df.tmp)=="origOverlap"] <- "overlap"
        df.tmp$distOv <- NULL

        df.out <- chutils::ch.rbind(df.out, df.tmp)
      }

      #combine the simulated data with the empirical data
      df.out1 <- df.out[,RWkeepColumns]
      df.out2 <- merge(df.out1, data, by.x = mergeByRWColumns, by.y=mergeByDataColumns)
    } else {

      if(!is.null(startBiasEffectCol)) {
        #if there is a start bias effect parameter, get the start bias effect conditions.  This should be either dummy coded or effect coded.
        sbConds <- unique(data[[startBiasEffectCol]])

        #now, for each start bias effect condition, run the simple random walk with the appropriate parameters and rbind them together
        df.out <- NULL
        for(sC in sbConds) {
          sIn <- startValue + (sC*startBiasEffect)
          df.tmp <- getMomentsOfRRWoverlap(overlapSeq, b, startValue = sIn, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops = loops, progress=progress)
          df.tmp[[startBiasEffectCol]] <- sC
          df.out <- chutils::ch.rbind(df.out, df.tmp)
        }

        #combine the simulated data with the empirical data
        df.out1 <- df.out[,RWkeepColumns]
        df.out2 <- merge(df.out1, data, by.x = mergeByRWColumns, by.y=mergeByDataColumns)
      } else {
        if(!is.null(valueChangeEffectCol)) {
          #if there is a start bias effect parameter, get the start bias effect conditions.  This should be either dummy coded or effect coded.
          vcConds <- unique(data[[valueChangeEffectCol]])

          #now, for each start bias effect condition, run the simple random walk with the appropriate parameters and rbind them together
          df.out <- NULL
          for(vC in vcConds) {
            overlapSeqIn <- overlapSeq + (vC * valueChangeEffect)
            df.tmp <- getMomentsOfRRWoverlap(overlapSeqIn, b, startValue = startValue, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops = loops, progress=progress)
            df.tmp[[valueChangeEffectCol]] <- vC

            #replace distorted overlaps (overlapSeq + vCp*valueChangeEffect) with the original overlaps, so it looks like real data
            df.ov <- data.frame (overlap = overlapSeqIn, origOverlap = overlapSeq)
            df.tmp <- merge (df.tmp, df.ov, by="overlap")
            colnames(df.tmp)[colnames(df.tmp)=="overlap"] <- "distOv"
            colnames(df.tmp)[colnames(df.tmp)=="origOverlap"] <- "overlap"
            df.tmp$distOv <- NULL

            df.out <- chutils::ch.rbind(df.out, df.tmp)

          }

          #combine the simulated data with the empirical data
          df.out1 <- df.out[,RWkeepColumns]
          df.out2 <- merge(df.out1, data, by.x = mergeByRWColumns, by.y=mergeByDataColumns)

        } else {
          df.out <- getMomentsOfRRWoverlap(overlapSeq, b, startValue = startValue, noiseSD = noiseSD, decayAsymptote = decayAsymptote, decayBeta = decayBeta, loops = loops, progress=progress)
          #combine the simulated data with the empirical data
          df.out1 <- df.out[,RWkeepColumns]
          df.out2 <- merge(df.out1, data, by.x = mergeByRWColumns, by.y=mergeByDataColumns)
        }
      }
    }

    #scale the RRW samples to fit the empirical RT data
    fit.lm <- NULL
    tryCatch ({
      fit.lm <- lm(df.out2[[dataRtCol]]~df.out2[[RwSamplesCol]])
    }, error = function(e) {})
    if(!is.null(fit.lm)) {
      df.out2$rtFit = df.out2[[RwSamplesCol]]*coef(fit.lm)[2]+coef(fit.lm)[1]
    } else {
      df.out2$rtFit = NA
    }

    return(df.out2)

}
