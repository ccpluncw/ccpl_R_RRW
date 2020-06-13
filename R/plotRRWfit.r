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
#' @param plotFilename A string that identifies the name of file (.pdf) in which the plot will be saved. The default is NULL, whereby the plot will not be saved.
#' @param twoPlotsPerPage A boolean that identifies whether to print two plots per page.  DEFAULT = TRUE.
#' @param yMinMixRT A vector of 2 numbers that identifies the c(min, max) for the y-axis of the RT graph. If not entered, then the function will calculate a pretty min and max. DEFAULT = NULL.
#' @return .
#' @keywords RRW random walk plot output fit
#' @export
#' @examples plotRRWFit (data, "rt", "pHit", "rtFit", "pHitFit", "correct", "overlap")

plotRRWFit <- function (data, dataRtCol = "rt", dataPhitCol = "pHit",  rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", condCol = NULL, plotFilename = NULL, twoPlotsPerPage = TRUE, yMinMixRT = NULL) {

    lwd <- 2

    if (twoPlotsPerPage == TRUE) {
      op <- par(mfrow=c(2,1),bty="n", font=1, family='serif', mar=c(2,5,2,5), oma=c(3,0,3,0), cex=1.25, las=1)
    } else {
      op <- par(mfrow=c(1,1), bg="white",  bty="n", font=2, family='serif', mar=c(5,6,4,10), las=1, cex=1)
      if (!is.null(plotFilename)) {
        pdf(plotFilename, width=10, height=8 )
      }
    }

    xlab = expression(paste("", Psi,"(value) Distributional overlap", sep=""))
    if(is.null(yMinMixRT)) {
      yMinMixRT <- chutils::ch.getPlotAxisMinMax(data[[dataRtCol]])
    }

    if(!is.null(condCol)) {

      conds <- unique(data[[condCol]])
      conds.n <- length(conds)
      for(cnd in 1:conds.n) {
        df.tmp <- data[ data[[condCol]]==conds[cnd] & data[[correctCol]] == TRUE & !is.na(data[[dataPhitCol]]), ]
        pty <- 16 + cnd - 1
        lty <- cnd
        if(cnd == 1) {
          plot(df.tmp[,overlapCol], df.tmp[,dataPhitCol], col=gray(0), pch=pty, ylim=c(0,1), xlab = xlab, ylab=NA)
          abline(a=0.5,b=0,col="grey", lwd=2)
        } else {
          points(df.tmp[,overlapCol], df.tmp[,dataPhitCol], col=gray(0), pch=pty)
        }
        lines(df.tmp[,overlapCol], df.tmp[, pHitFitCol], col=gray(0), lty=lty, lwd = lwd)
      }

      for(cnd in 1:conds.n) {
        df.tmp <- data[ data[[condCol]]==conds[cnd] & !is.na(data[[dataRtCol]]), ]
        pty <- 16 + cnd - 1
        lty <- cnd
        if(cnd == 1) {
          plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=pty, ylim=yMinMixRT, xlab = xlab, ylab=NA)
        } else {
          points(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=pty)
        }
        points(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE,dataRtCol], col=gray(.6), pch=pty)
        lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, rtFitCol], col=gray(0), lty=lty, lwd = lwd)
        lines(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE, rtFitCol], col=gray(.6), lty=lty, lwd = lwd)
      }

    } else {

          df.tmp <- data[ !is.na(data[[dataPhitCol]]), ]

          plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataPhitCol], col=gray(0), pch=16, ylim=c(0,1), xlab = xlab, ylab=NA)
          abline(a=0.5,b=0,col="grey", lwd=2)
          lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, pHitFitCol], col=gray(0), lty='solid', lwd = lwd)

          df.tmp <- data[ !is.na(data[[dataRtCol]]), ]

          plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=16, ylim=yMinMixRT, xlab = xlab, ylab=NA)
          points(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE,dataRtCol], col=gray(.6), pch=16)
          lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, rtFitCol], col=gray(0), lty='solid', lwd = lwd)
          lines(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE, rtFitCol], col=gray(.6), lty='solid', lwd = lwd)
    }

    if (!is.null(plotFilename) & twoPlotsPerPage == TRUE) {
			  dev.copy(pdf, plotFilename, width=6, height=9)
				dev.off();
    }
    par(op)
}
