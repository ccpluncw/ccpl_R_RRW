#' This function plots the fitted RRW simulation as a function of the empirical data.
#'
#' Function plots the fitted RRW simulation as a function of the empirical data.  It is generally run after the getRRWfit().
#'
#'
#'
#' @param data This is a dataframe that must contain the following columns: overlap; RT (often a median); the proportion correct/incorrect; whether or not the row specifies a correct or incorrect trial. The dataset can also contain a column specifying a condition that will influence either the startpoint or the value of the trials.
#' @param dataRtCol A string that identifies the name of the column in data that contains the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param dataPhitCol A string that identifies the name of the column in data that contains the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param rtFitCol A string that identifies the name of the column in data that contains the RRW fit of the RTs for the specific overlap/correct/condition combination. Default is "rt"
#' @param pHitFitCol A string that identifies the name of the column in data that contains the RRW fit of the proportion of trials that are either correct or incorrect for the specific overlap/correct/condition combination. Default is "pHit"
#' @param correctCol A string that identifies the name of the column in data that identifies whether the trials were correct (TRUE) or incorrect (FALSE). The default is "correct"
#' @param overlapCol A string that identifies the name of the column in data that contains the distributional overlaps for each row. The default is "overlap"
#' @param startBiasEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the position of the start point.  This column should be effect coded, because the values contained in this column will be multiplied by the value of startBiasEffect (and then added to s) to get the start point of the RRW.
#' @param valueChangeEffectCol A string that identifies the name of the column in data that identifies the conditions that will affect the change of overlap (value change).  This column should be effect coded, because the values contained in this column will be multiplied by the value of valueChangeEffect (and then added to overlap) to get the overlap input in the RRW. This change is assumed to be constant across the entire overlap sequence (from 0-1).
#' @param plotFilename A string that identifies the name of file (.pdf) in which the plot will be saved. The default is NULL, whereby the plot will not be saved.
#' @param twoPlotsPerPage A boolean that identifies whether to print two plots per page.  DEFAULT = TRUE.
#'
#' @return .
#'
#' @keywords RRW random walk plot output fit
#' @export
#' @examples plotRRWFit (data, "rt", "pHit", "rtFit", "pHitFit", "correct", "overlap")

plotRRWFit <- function (data, dataRtCol, dataPhitCol,  rtFitCol, pHitFitCol, correctCol, overlapCol, startBiasEffectCol = NULL, valueChangeEffectCol = NULL,  plotFilename, twoPlotsPerPage = TRUE) {

    if (twoPlotsPerPage == TRUE) {
      op <- par(mfrow=c(2,1),bty="n", font=1, family='serif', mar=c(2,5,2,5), oma=c(3,0,3,0), cex=1.25, las=1)
    } else {
      op <- par(mfrow=c(1,1), bg="white",  bty="n", font=2, family='serif', mar=c(5,6,4,10), las=1, cex=1)
      if (!is.null(plotFilename)) {
        pdf(plotFilename, width=10, height=8 )
      }
    }

    xlab = expression(paste("", Psi,"(value) Distributional overlap", sep=""))
    yMinMax <- chutils::ch.getPlotAxisMinMax(data[[dataRtCol]])

    if(!is.null(startBiasEffectCol)) {
      sbConds <- unique(data[[startBiasEffectCol]])

      for(sC in 1:length(sbConds)) {
        df.tmp <- data[data[[startBiasEffectCol]]==sbConds[sC],]
        pty <- 16 + sC - 1
        lty <- sC
        if(sC == 1) {
          plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataPhitCol], col=gray(0), pch=pty, ylim=c(0,1), xlab = xlab, ylab=NA)
          abline(a=0.5,b=0,col="grey", lwd=2)
        } else {
          points(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataPhitCol], col=gray(0), pch=pty)
        }
        lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, pHitFitCol], col=gray(0), lty=lty)
      }

      for(sC in 1:length(sbConds)) {
        df.tmp <- na.omit(data[data[[startBiasEffectCol]]==sbConds[sC],])
        pty <- 16 + sC - 1
        lty <- sC
        if(sC == 1) {
          plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=pty, ylim=yMinMax, xlab = xlab, ylab=NA)
        } else {
          points(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=pty)
        }
        points(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE,dataRtCol], col=gray(.6), pch=pty)
        lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, rtFitCol], col=gray(0), lty=lty)
        lines(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE, rtFitCol], col=gray(.6), lty=lty)
      }

    } else {

      if(!is.null(valueChangeEffectCol)) {
        vcConds <- unique(data[[valueChangeEffectCol]])

        for(vC in 1:length(vcConds)) {
          df.tmp <- data[data[[valueChangeEffectCol]]==vcConds[vC],]
          pty <- 16 + vC - 1
          lty <- vC
          if(vC == 1) {
            plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataPhitCol], col=gray(0), pch=pty, ylim=c(0,1), xlab = xlab, ylab=NA)
            abline(a=0.5,b=0,col="grey", lwd=2)
          } else {
            points(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataPhitCol], col=gray(0), pch=pty)
          }
          lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, pHitFitCol], col=gray(0), lty=lty)
        }

        for(vC in 1:length(vcConds)) {
          df.tmp <- na.omit(data[data[[valueChangeEffectCol]]==vcConds[vC],])
          pty <- 16 + vC - 1
          lty <- vC
          if(vC == 1) {
            plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=pty, ylim=yMinMax, xlab = xlab, ylab=NA)
          } else {
            points(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=pty)
          }
          points(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE,dataRtCol], col=gray(.6), pch=pty)
          lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, rtFitCol], col=gray(0), lty=lty)
          lines(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE, rtFitCol], col=gray(.6), lty=lty)
        }

      } else {
        df.tmp <- na.omit(data)

        plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataPhitCol], col=gray(0), pch=16, ylim=c(0,1), xlab = xlab, ylab=NA)
        abline(a=0.5,b=0,col="grey", lwd=2)
        lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, pHitFitCol], col=gray(0), lty='solid')


        plot(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE,dataRtCol], col=gray(0), pch=16, ylim=yMinMax, xlab = xlab, ylab=NA)
        points(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE,dataRtCol], col=gray(.6), pch=16)
        lines(df.tmp[df.tmp[[correctCol]] == TRUE,overlapCol], df.tmp[df.tmp[[correctCol]] == TRUE, rtFitCol], col=gray(0), lty='solid')
        lines(df.tmp[df.tmp[[correctCol]] == FALSE,overlapCol], df.tmp[df.tmp[[correctCol]] == FALSE, rtFitCol], col=gray(.6), lty='solid')
      }
    }

    if (!is.null(plotFilename) & twoPlotsPerPage == TRUE) {
			  dev.copy(pdf, plotFilename, width=6, height=9)
				dev.off();
    }
    par(op)
}
