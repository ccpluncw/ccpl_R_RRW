#' This function plots the rrw data fits from the output of getMeanRRWfit() and the rrwModelList.
#'
#' Function that plots the rrw data fits from the output of getMeanRRWfit() and the rrwModelList.
#' @param df.RRWout This is a dataframe containing the output from getMeanRRWfit() or rrwRunWithFixedParameters () or rrwRunSmartGridSearch (). It contains the rrwData and the simulation fits.
#' @param rrwModelList A list that specifies the rrw model.  Build the rrwModelList useing rrwAddEffectToRRWModel
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param rtFitCol A string that identifies the name of the column in data that contains the RRW fit of the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param pHitFitCol A string that identifies the name of the column in data that contains the RRW fit of the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pCross"
#' @param correctCol a boolean that specifies the name of the column that specifies if response was correct or incorrect.  It must be coded as a boolean, whereby correct == TRUE, and incorrect == FALSE. The default is "correct"
#' @param overlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param fileTag A string that is appended to the name of files to identify the analysis and experiment. The default is NULL, whereby the filetag will just be based on a timestamp.
#' @param numSimsToPlot A number specifying how many simulation runs are in the dataset and should be plotted. Default is 10.
#' @param maxIntensityChanges the maximum number of distinguishable intensity changes. DEFAULT = 8.
#' @param maxHueChanges the maximum number of distinguishable hue changes. DEFAULT = 10.
#''
#' @return Several files containing plots are saved.  The plots are duplicated: one version with each simulation plus the average and one with just the average.
#' @keywords RRW plots
#' @export
#' @examples rrwPlotSGSoutput(df.fitted, "rt", "pHit", "rtFit", "pCross", "correct", "overlap", fileTag = NULL, numSimsToPlot = 40)

rrwPlotSGSoutput2 <- function (df.RRWout, rrwModelList, dataRtCol = 'rt', dataPhitCol = 'pHit', rtFitCol = 'rtFit', pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap",fileTag = NULL, numSimsToPlot = 10, maxIntensityChanges = 8, maxHueChanges = 10) {

  #if fileTag does not exist, create one based on the timestamp
  if(is.null(fileTag)) {
    fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), "RRW", sep="_")
  }

  #create the plotname
  plotFilename <- paste(fileTag,"FitResults.pdf")
  #get the xy axis bounds
  yMinMixRT <- chutils::ch.getPlotAxisMinMax(df.RRWout[[dataRtCol]])

  #add the grouping columns specified in the rrwModelList
  grpCols <- NULL
  #try as a nested list first
  tryCatch ({
    grpCols <- unique(as.vector(unlist(lapply(unlist(rrwModelList, recursive=FALSE), `[`, "GroupByVariables"))))
    }, error = function(e) {}
  )
  #if that fails, try as an unnested list
  if(is.null(grpCols)) {
    tryCatch ({
      grpCols <- unique(as.vector(unlist(lapply(rrwModelList, `[`, "GroupByVariables"))))
      }, error = function(e) {}
    )
  }

  #create a condition variable that codes all the conditions in a single variable
  #if there are grouping variables
  if(!is.null(grpCols)) {
    #how many grouping variables
    nGrpsCols <- length(grpCols)
    if( nGrpsCols > 3) {
      #for each group
      for (i in grpCols) {
        #if it is the first group
        if(i == grpCols[1]) {
          #initiate the variable
          df.RRWout$cond <- df.RRWout[[i]]
        } else {
          #otherise, just add values to the variable
          df.RRWout$cond <- paste(df.RRWout$cond, df.RRWout[[i]], sep="-")
        }
      }
      df.legend <- ch.getPlotLegendValsFromData(df.RRWout,"cond", maxIntensityChanges = maxIntensityChanges, maxHueChanges = maxHueChanges)
      grpCols <- "cond"
    } else {
      df.legend <- ch.getPlotLegendValsFromData(df.RRWout,grpCols, maxIntensityChanges = maxIntensityChanges, maxHueChanges = maxHueChanges)
    }
  }
  #for each grouping variable, plot the data with all the simulations and without all the simulations (only the mean simulation).
  if(!is.null(grpCols)) {
    nGrpsCols <- length(grpCols)
    for (i in grpCols) {
      tmpConds <- unique(df.RRWout[[i]])
      for(j in tmpConds) {
        df.legend.tmp <- df.legend[df.legend[[i]] == j,]
        #plot with individual simulations
        plotRRWFit3(df.RRWout[df.RRWout[[i]] == j,], dataRtCol, dataPhitCol,  rtFitCol, pHitFitCol, correctCol, overlapCol, condCol =grpCols,  plotFilename = paste("withSims-",i,"-",j, plotFilename, sep=""), yMinMixRT=yMinMixRT, numSimsToPlot=numSimsToPlot, df.legend = df.legend.tmp)
        #plot without individual simulations
        plotRRWFit3(df.RRWout[df.RRWout[[i]] == j,], dataRtCol, dataPhitCol,  rtFitCol, pHitFitCol, correctCol, overlapCol, condCol =grpCols,  plotFilename = paste(i,"-",j, plotFilename, sep=""), yMinMixRT=yMinMixRT, numSimsToPlot=0, df.legend = df.legend.tmp)
      }
    }
  }

  #for the condition variable, plot the data with all the simulations and without all the simulations (only the mean simulation).
    #plot with individual simulations
    plotRRWFit3(df.RRWout, dataRtCol, dataPhitCol,  rtFitCol, pHitFitCol, correctCol, overlapCol, condCol =grpCols,  plotFilename = paste("withSim-simple", plotFilename), yMinMixRT=yMinMixRT, numSimsToPlot=numSimsToPlot, df.legend = df.legend)
    #plot without individual simulations
    plotRRWFit3(df.RRWout, dataRtCol, dataPhitCol,  rtFitCol, pHitFitCol, correctCol, overlapCol, condCol =grpCols,  plotFilename = paste("simple", plotFilename), yMinMixRT=yMinMixRT, numSimsToPlot=0, df.legend = df.legend)
}
