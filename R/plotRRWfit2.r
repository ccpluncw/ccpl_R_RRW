#' This function plots the fitted RRW simulation as a function of the empirical data.
#'
#' Function plots the fitted RRW simulation as a function of the empirical data.  It is generally run after the getRRWfit().
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contains columns that effect code the influence of different parameters.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param rtFitCol A string that identifies the name of the column in data that contains the RRW fit of the RTs for the specific overlap/correct/condition combination. Default is "rtFit"
#' @param pHitFitCol A string that identifies the name of the column in data that contains the RRW fit of the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pCross"
#' @param correctCol A string that identifies the name of the column in data that identifies whether the trials were correct (TRUE) or incorrect (FALSE). The default is "correct"
#' @param overlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param condCol A string that identifies the name of the column in data that identifies the conditions that will be plotted as separate lines. DEFAULT = NULL
#' @param numSimsToPlot A number specifying how many simulation runs are in the dataset and should be plotted. Default is 0, indicating that no simulations should be plotted.
#' @param plotFilename A string that identifies the name of file (.pdf) in which the plot will be saved. The default is NULL, whereby the plot will not be saved.
#' @param multiplePlotsPerPage A boolean that identifies whether to print multiple plots per page.  DEFAULT = TRUE.
#' @param yMinMixRT A vector of 2 numbers that identifies the c(min, max) for the y-axis of the RT graph. If not entered, then the function will calculate a pretty min and max. DEFAULT = NULL.
#' @param xMinMixRT A vector of 2 numbers that identifies the c(min, max) for the x-axis (overlap). DEFAULT = c(0,1).
#' @param combineRThvoRTlvo A boolean that specifies whether to print the RT_HVO and RT_LVO on the same graph (TRUE) or separate graphs (FALSE).  DEFAULT = FALSE.
#' @param maxIntensityChanges the maximum number of distinguishable intensity changes. DEFAULT = 8.
#' @param maxHueChanges the maximum number of distinguishable hue changes. DEFAULT = 10.
#' @return .
#' @keywords RRW random walk plot output fit
#' @export
#' @examples plotRRWFit2 (data, "rt", "pHit", "rtFit", "pHitFit", "correct", "overlap")

plotRRWFit2 <- function (data, dataRtCol = "rt", dataPhitCol = "pHit",  rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", condCol = NULL, numSimsToPlot = 0, plotFilename = NULL, multiplePlotsPerPage = TRUE, yMinMixRT = NULL, xMinMax = NULL, combineRThvoRTlvo = FALSE, maxIntensityChanges = 8, maxHueChanges = 10) {


  #set up plot parameters
    if (multiplePlotsPerPage == TRUE) {
      op <- par(mfrow=c(2,2),bty="n", font=1, family='serif', mar=c(2,2,2,2), oma=c(2,2,2,2), cex=1.25, las=1)
    } else {
      op <- par(mfrow=c(1,1), bg="white",  bty="n", font=2, family='serif', mar=c(5,6,4,10), las=1, cex=1)
      if (!is.null(plotFilename)) {
        pdf(plotFilename, width=10, height=8 )
      }
    }
    lwd <- 2
    pchTmp <- 21
    simGrey <- 0.95
    pointCEXsize <- 0.7


  #set up plot axes
    if(is.null(yMinMixRT)) {
      yMinMixRT <- chutils::ch.getPlotAxisMinMax(data[[dataRtCol]])
    }
    if(is.null(xMinMax)) {
      minX <- ifelse(min(data[[overlapCol]]) < 0, min(data[[overlapCol]]), 0)
      maxX <- ifelse(max(data[[overlapCol]]) > 1, max(data[[overlapCol]]), 1.0)
      xLims <- c(minX, maxX)
    } else {
      xLims <- xMinMax
    }


  #get conditions and number of conditions
    if(!is.null(condCol)) {
      conds <- unique(data[[condCol]])
      conds.n <- length(conds)
      #create plot info and legend info
      df.grpIdx <- data.frame(cond = conds, idx = seq(1,conds.n, 1))
      df.legend <- chutils::ch.getPlotLegendVals(df.grpIdx, maxIntensityChanges = maxIntensityChanges, maxHueChanges = maxHueChanges)
    } else {
      conds.n <- 1
      df.grpIdx <- data.frame(cond = "all", idx = c(1))
      df.legend <- chutils::ch.getPlotLegendVals(df.grpIdx, maxIntensityChanges = maxIntensityChanges, maxHueChanges = maxHueChanges)
    }

    #create empty plot for p(HVO)
    plot(NA, ylim=c(0,1), xlim = xLims, xlab = NULL, ylab= NA)

    abline(a=0.5,b=0,col="grey", lwd=2)

    if(numSimsToPlot > 0) {
      for(cnd in 1:conds.n) {
        if(!is.null(condCol)) {
          df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == TRUE & !is.na(data[[dataPhitCol]]), ]
        } else {
          df.tmp <- data[ data[[correctCol]] == TRUE & !is.na(data[[dataPhitCol]]), ]
        }

        plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
        lty <- df.legend[cnd,"lty"]

        #add background Simulation lines
          for(i in 1: numSimsToPlot) {
            lines(df.tmp[,overlapCol], df.tmp[, paste(pHitFitCol,i,sep="")], col=gray(simGrey), lty="solid", lwd = 1)
          }
        }
      }

      for(cnd in 1:conds.n) {
        if(!is.null(condCol)) {
          df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == TRUE & !is.na(data[[dataPhitCol]]), ]
        } else {
          df.tmp <- data[ data[[correctCol]] == TRUE & !is.na(data[[dataPhitCol]]), ]
        }

        plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
        lty <- df.legend[cnd,"lty"]

        #add data points and final fit
        points(df.tmp[[overlapCol]], df.tmp[[dataPhitCol]], col=plotCol, pch=pchTmp, bg = plotCol, cex=pointCEXsize)
        lines(df.tmp[[overlapCol]], df.tmp[[pHitFitCol]], col=plotCol, lty=lty, lwd = lwd)
      }

    #create an empty plot for the legend
    plot.new()
    chutils::ch.addLegend(df.legend, "cond", placement="center", cexLegend = 0.75, includeTitle = F)

    #create empty plot for RT_HVO
    plot(NA, ylim=yMinMixRT, xlim = xLims, xlab = NULL, ylab=NA)

      if(numSimsToPlot > 0) {
        for(cnd in 1:conds.n) {
          if(!is.null(condCol)) {
            df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == TRUE & !is.na(data[[dataRtCol]]), ]
          } else {
            df.tmp <- data[ data[[correctCol]] == TRUE & !is.na(data[[dataRtCol]]), ]
          }

          plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
          lty <- df.legend[cnd,"lty"]

          #add background Simulation lines
          for(i in 1: numSimsToPlot) {
            rtColTmp <- paste(rtFitCol,i,sep="")
            lines(df.tmp[[overlapCol]], df.tmp[[rtColTmp]], col=gray(simGrey), lty="solid", lwd = 1)
          }
        }
      }

      #if combining HVO and LVO, then draw LVO simulations before you draw the HVO datapoints
      if(combineRThvoRTlvo == TRUE) {
        if(numSimsToPlot > 0) {
          for(cnd in 1:conds.n) {
            if(!is.null(condCol)) {
              df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == FALSE & !is.na(data[[dataRtCol]]), ]
            } else {
              df.tmp <- data[ data[[correctCol]] == FALSE & !is.na(data[[dataRtCol]]), ]
            }

            plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
            lty <- df.legend[cnd,"lty"]

            #add background Simulation lines
            for(i in 1: numSimsToPlot) {
              rtColTmp <- paste(rtFitCol,i,sep="")
              lines(df.tmp[[overlapCol]], df.tmp[[rtColTmp]], col=gray(simGrey), lty="solid", lwd = 1)
            }
          }
        }
      }

      for(cnd in 1:conds.n) {
        if(!is.null(condCol)) {
          df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == TRUE & !is.na(data[[dataRtCol]]), ]
        } else {
          df.tmp <- data[ data[[correctCol]] == TRUE & !is.na(data[[dataRtCol]]), ]
        }

        plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
        lty <- df.legend[cnd,"lty"]

        #add data points and final fit
        points(df.tmp[[overlapCol]], df.tmp[[dataRtCol]], col=plotCol, pch=pchTmp, bg = plotCol, cex=pointCEXsize)
        lines(df.tmp[[overlapCol]], df.tmp[[rtFitCol]], col=plotCol, lty=lty, lwd = lwd)
      }

    #if you are not combining HVO and LVO, then create empty plot for RT_LVO
    if(combineRThvoRTlvo == FALSE) {
      plot(NA, ylim=yMinMixRT, xlim = xLims, xlab = NULL, ylab=NA)

      #draw the LVO simulations now
      if(numSimsToPlot > 0) {
        for(cnd in 1:conds.n) {
          if(!is.null(condCol)) {
            df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == FALSE & !is.na(data[[dataRtCol]]), ]
          } else {
            df.tmp <- data[ data[[correctCol]] == FALSE & !is.na(data[[dataRtCol]]), ]
          }

          plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
          lty <- df.legend[cnd,"lty"]

          #add background Simulation lines
          for(i in 1: numSimsToPlot) {
            rtColTmp <- paste(rtFitCol,i,sep="")
            lines(df.tmp[[overlapCol]], df.tmp[[rtColTmp]], col=gray(simGrey), lty="solid", lwd = 1)
          }
        }
      }
    } else {
      #if you are combining the HVO and LVO, then change the LVO symbol to an X
      pchTmp <- 4
      pointCEXsize <- 0.5
    }

      for(cnd in 1:conds.n) {
        if(!is.null(condCol)) {
          df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == FALSE & !is.na(data[[dataRtCol]]), ]
        } else {
          df.tmp <- data[ data[[correctCol]] == FALSE & !is.na(data[[dataRtCol]]), ]
        }

        plotCol <- hsv(df.legend[cnd,"h"], df.legend[cnd,"s"], df.legend[cnd,"v"])
        lty <- df.legend[cnd,"lty"]

        #add data points and final fit
        points(df.tmp[[overlapCol]], df.tmp[[dataRtCol]], col=plotCol, pch=pchTmp, bg = plotCol, cex=pointCEXsize)
        lines(df.tmp[[overlapCol]], df.tmp[[rtFitCol]], col=plotCol, lty=lty, lwd = lwd)
      }

    if (!is.null(plotFilename) & multiplePlotsPerPage == TRUE) {
			  dev.copy(pdf, plotFilename, width=10, height=8)
				dev.off();
    }
    par(op)
}
