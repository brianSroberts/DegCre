#' Plot DegCre Bin Algorithm Statistics
#'
#' Plots the DegCre distance bin optimization statistic against different bin
#' sizes, highlighting the optimal bin size.
#'
#' @param degCreResList A list of DegCre results.
#'
#' @details
#' This function takes a DegCre results list and plots the bin heuristic
#' statistics against different bin sizes. It also highlights the optimal bin
#' size chosen based on the analysis.
#' The y-axis of the plot is the median KS statistic of all bins versus the
#' global CRE p-value distribution.
#'
#' @seealso \link{distBinHeuristic} for calculating the DEG-CRE bin heuristic.
#'
#' @return Invisibly, the picked optimal bin size.
#'
#' @examples
#' #Load example data.
#' data(DexNR3C1)
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
#'                                    DegP=DexNR3C1$DegGR$pVal,
#'                                    DegLfc=DexNR3C1$DegGR$logFC,
#'                                    CreGR=DexNR3C1$CreGR,
#'                                    CreP=DexNR3C1$CreGR$pVal,
#'                                    CreLfc=DexNR3C1$CreGR$logFC)
#'
#' #Plot distance bin median KS statistic curve.
#'
#'plotDegCreBinHeuristic(degCreResList=degCreResListDexNR3C1)
#'
#' @author Brian S. Roberts
#'
#' @export
plotDegCreBinHeuristic <- function(degCreResList) {
  binSizeKSstatMat <- degCreResList$binHeurOutputs$crePKsMat
  pickedBinSize <- degCreResList$binHeurOutputs$pickedBinSize
  maskPicked <- which(binSizeKSstatMat[, 1] == pickedBinSize)

  plot(x = binSizeKSstatMat[, 1],
       y = binSizeKSstatMat[, 2],
       xlab = "Bin Size (number of hits)",
       ylab = "Median bin KS stat",
       type = 'b',
       pch = 16,
       log = 'x',
       cex = 0.8)

  points(x = binSizeKSstatMat[maskPicked, 1],
         y = binSizeKSstatMat[maskPicked, 2],
         pch = 16,
         col = "red",
         cex = 0.8)

  yRange <- max(binSizeKSstatMat[, 2]) - min(binSizeKSstatMat[, 2])
  yText <- min(binSizeKSstatMat[, 2]) + 0.85 * yRange

  xRange <- max(binSizeKSstatMat[, 1]) - min(binSizeKSstatMat[, 1])
  xText <- min(binSizeKSstatMat[, 1]) + 0.25 * xRange

  par(xpd = NA)
  text(x = xText, y = yText,
    label = paste("optimal bin size", pickedBinSize, sep = "\n"),
    col = "red",cex=0.6)
  invisible(binSizeKSstatMat[maskPicked,1])
}

#' Plot DegCre Association Probability vs. Binned Genomic Distance
#'
#' Plots the DegCre association probability against binned genomic distance and
#' highlights the quantile range.
#'
#' @param degCreResList A list of DegCre results.
#' @param assocProbFDRThresh Numeric value from 0 to 1 specifying the FDR
#' threshold for association probability. (Default: \code{0.05})
#' @param plotQRange Numeric vector of quantile range for plotting
#' (e.g., \code{c(0.25, 0.75)} for interquartile range).
#' (Default: \code{c(0.25, 0.75)})
#' @param hiYLim Numeric value specifying the upper limit of the y-axis.
#' (Default: \code{NULL})
#' @param loYLim Numeric value specifying the lower limit of the y-axis.
#' (Default: \code{NULL})
#' @param qRangeFillColor Color for filling the quantile range polygon.
#' (Default: \code{#88CCEE})
#' @param nullLineColor Color for the null association probability line.
#' (Default: \code{#CC6677})
#' @param plotLabCex Numeric value specifying the character expansion factor
#' for labels. (Default: 0.8)
#' @param plotAxisCex Numeric value specifying the character expansion factor
#' for axis labels. (Default: 0.5)
#'
#' @details
#' This function takes the results of the DegCre analysis, including genomic
#' distances and association probabilities, and creates a plot of association
#' probabilities against binned genomic distances. It highlights the quantile
#' range (e.g., interquartile range) and includes a line for null association
#' probabilities.
#' The top panel shows the number of associations passing
#' \code{assocProbFDRThresh}. The bottom panel shows the median FDR-passing
#' association probability as a black line, with the specified quantile range
#' (defaults to interquartile) plotted as \code{qRangeFillColor} region.
#' The \code{nullLineColor} colored line is the null association probability,
#' that is the association probability for a bin with uniform CRE p-values.
#
#' @return Invisibly, a matrix with these columns:
#' \describe{
#'   \item{binMidDist}{Numeric value of the midpoint distance of the bin
#'   (TSS to CRE) in kb.}
#'   \item{q_<\code{plotQRange[1] x 100}>}{Numeric value of lower bound of the
#'   highlight region.}
#'   \item{q_50}{Numeric value of plotted line.}
#'   \item{q_<\code{plotQRange[2] x 100}>}{Numeric value of upper bound of the
#'   highlight region.}
#'   \item{nullAssocProb}{Numeric of null association probability of the bin.}
#' }
#'
#' @examples
#' #Load example data.
#' data(DexNR3C1)
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
#'                                    DegP=DexNR3C1$DegGR$pVal,
#'                                    DegLfc=DexNR3C1$DegGR$logFC,
#'                                    CreGR=DexNR3C1$CreGR,
#'                                    CreP=DexNR3C1$CreGR$pVal,
#'                                    CreLfc=DexNR3C1$CreGR$logFC)
#'
#' #Plot association probability versus binned genomic distance.
#'
#' outProbVsDistMat <-
#'  plotDegCreAssocProbVsDist(degCreResList=degCreResListDexNR3C1)
#'
#' @author Brian S. Roberts
#'
#' @export
plotDegCreAssocProbVsDist <- function(degCreResList,
                                     assocProbFDRThresh = 0.05,
                                     plotQRange = c(0.25, 0.75),
                                     hiYLim = NULL,
                                     loYLim = NULL,
                                     qRangeFillColor = "#88CCEE",
                                     nullLineColor = "#CC6677",
                                     plotLabCex = 0.6,
                                     plotAxisCex = 0.5) {

    hitsDegCre <- degCreResList$degCreHits
    dists <- mcols(hitsDegCre)$assocDist

    dists <- dists/1000

    assocProbs <- mcols(hitsDegCre)$assocProb
    assocProbFDRs <- mcols(hitsDegCre)$assocProbFDR
    maxBinDists <- mcols(hitsDegCre)$binAssocDist

    maxBinDists <- maxBinDists/1000

    uniqMaxBinDists <- sort(unique(maxBinDists))

    maskPassFDR <- which(assocProbFDRs<=assocProbFDRThresh)

    plotDists <- dists

    inParMai <- par("mai")
    inParBty <- par("bty")

    par(mai=c(0,0,0,0))

    layX <- layout(matrix(seq_len(2),ncol=1),heights=c(50,100))

    firstPlotMai <- inParMai
    firstPlotMai[1] <- 0

    par(mai=firstPlotMai)

    firstPlotBty <- "n"
    par(bty=firstPlotBty)

    #plot text of how many pass threshold in each bin
    numPassPerBin <- unlist(lapply(uniqMaxBinDists,function(binDistX){
        maskInBinX <- which(maxBinDists==binDistX)
        return(length(intersect(maskInBinX,maskPassFDR)))
    }))

    #get bin midpoints
    binMids <- c(0,uniqMaxBinDists[seq_len((length(uniqMaxBinDists)-1))]) +
        diff(c(0,uniqMaxBinDists))/2

    xLimits <- c(0,max(plotDists))

    nPassYlab <- "Num.\nPer Bin"


    par(xpd=NA)

    plot(x=binMids,y=numPassPerBin,type="l",lwd=1,col="black",
        ylab=nPassYlab,xlab="",xaxt='n',cex.lab=plotLabCex,
        ylim=c(0,1.02*max(numPassPerBin)),xlim=xLimits,cex.axis=plotAxisCex)

    #get median and quantile range  for all pass data and plot polygons

    quantsMat <- matrix(ncol=4)

    for(i in seq_along(uniqMaxBinDists)){
        binDistX <- uniqMaxBinDists[i]
        maskDistBinX <- which(maxBinDists==binDistX)

        maskPassDistBinX <- intersect(maskDistBinX,maskPassFDR)

        outQsX <- quantile(assocProbs[maskPassDistBinX],
            probs=c(plotQRange[1],0.5,plotQRange[2]))

        outRowX <- c(binDistX,outQsX)

        quantsMat <- rbind(quantsMat,outRowX)
    }

    quantsMat <- quantsMat[2:nrow(quantsMat),]

    lowerQName <- paste("q",plotQRange[1]*100,sep="_")
    upperQName <- paste("q",plotQRange[2]*100,sep="_")
    colnames(quantsMat) <- c("binMidDist",lowerQName,"q50",upperQName)

    #get null (distance based only) association probability
    nullAssocMat <- getDistBinNullAssocProb(degCreResList)
    nullAssocMat[,1] <- nullAssocMat[,1]/1000

    #replace NA rows in quantsMat with corresponding nullAssocMat val

    maskNARowQuants <- which(is.na(quantsMat[,2]))

    quantsMat[maskNARowQuants,2:4] <- nullAssocMat[maskNARowQuants,c(2,2,2)]

    allAssocs <- c(as.numeric(quantsMat[,2:4]),as.numeric(nullAssocMat[,2]))

    maxy <- max(allAssocs)
    if(is.null(hiYLim)){
        hiYLim <- 1.02*maxy
    }

    if(is.null(loYLim)){
        loYLim <- 0.98*min(allAssocs)
    }

    yLimits <- c(loYLim,hiYLim)

    secondPlotMai <- inParMai
    secondPlotMai[3] <- 0.15

    par(mai=secondPlotMai)
    par(bty=inParBty)

    xlabel <- "Bin Upper Dist. (Kb)"

    par(xpd=FALSE)
    par(cex.axis=plotAxisCex)
    par(cex.lab=plotLabCex)

    #initialize plot area
    plot(x=quantsMat[,1],y=quantsMat[,4],type='n',ylim=yLimits,
        ylab="Assoc.\nProbs.",xlab=xlabel,xaxt='n')

     tempMgp <- par('mgp')
     tempMgp[2] <- 0.5*tempMgp[2]
     tempMgp[1] <- 0.5*tempMgp[1]
     par(mgp=tempMgp)

     axis(side=1,cex.axis=plotAxisCex,cex.lab=plotLabCex)

    #make polygon of quantile range
    polyXs <- c(quantsMat[,1],quantsMat[(nrow(quantsMat):1),1])
    polyYs <- c(quantsMat[,2],quantsMat[(nrow(quantsMat):1),4])

    qRangeFill <- changeColorAlpha(qRangeFillColor,newAlpha=150)
    qRangeOutline <- qRangeFillColor

    polygon(x=polyXs,y=polyYs,col=qRangeFill,border=qRangeOutline)

    #plot median line
    lines(x=quantsMat[,1],y=quantsMat[,3],col="black")

    #plot null association prob
    lines(x=nullAssocMat[,1],y=nullAssocMat[,2],
        lwd=1.2,col=nullLineColor)

    outMat <- cbind(quantsMat,nullAssocMat[,2])
    colnames(outMat)[5] <- "nullAssocProb"
    rownames(outMat) <- NULL

    invisible(outMat)
}

#' Calculate Null Association Probability for Each Distance Bin
#'
#' Calculates the null association probability for each distance bin in the
#' DegCre analysis.
#'
#' @param degCreResList A list of DegCre results.
#'
#' @details
#' This function takes the results of the DegCre analysis and computes the null
#' association probability for each unique distance bin. The null association
#' probability represents the expected proportion of differentially expressed
#' genes (DEGs) in each distance bin under the null hypothesis.
#'
#' @return A matrix with these columns:
#' \describe{
#'   \item{binAssocDist}{Numeric value representing the distance bin
#'   (TSS to CRE) in base pairs.}
#'   \item{nullAssocProb}{Numeric value representing the null association
#'   probability of the bin.}
#' }
#'
#' @examples
#' #' #Load example data.
#' data(DexNR3C1)
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
#'                                    DegP=DexNR3C1$DegGR$pVal,
#'                                    DegLfc=DexNR3C1$DegGR$logFC,
#'                                    CreGR=DexNR3C1$CreGR,
#'                                    CreP=DexNR3C1$CreGR$pVal,
#'                                    CreLfc=DexNR3C1$CreGR$logFC)
#'
#' # Calculate null association probabilities.
#' outNullMat <- getDistBinNullAssocProb(degCreResList = degCreResListDexNR3C1)
#'
#' @author Brian S. Roberts
#'
#' @export
getDistBinNullAssocProb <- function(degCreResList) {
    hitsDegCre <- degCreResList$degCreHits

    binAssocDists <- mcols(hitsDegCre)$binAssocDist
    uniqBinAssocDists <- sort(unique(binAssocDists))

    alphaDEG <- degCreResList$alphaVal

    degPadjs <- mcols(hitsDegCre)$DegPadj

    binNullAssoc <- unlist(lapply(uniqBinAssocDists,function(uniqDistBinX){
        maskDistBinX <- which(binAssocDists==uniqDistBinX)
        nDegs <-
            length(which(degPadjs[maskDistBinX] <= alphaDEG)) * (1-alphaDEG)
        nullAssoc <- nDegs/length(maskDistBinX)
        return(nullAssoc)
    }))

    outMat <- cbind(uniqBinAssocDists,binNullAssoc)
    colnames(outMat) <- c("binAssocDist","nullAssocProb")
    return(outMat)
}

#' Change Color Transparency
#'
#' Changes the transparency (alpha channel) of a color or vector of colors.
#'
#' @param colorVec Character vector of hexadecimal or named colors to be
#' modified.
#' @param newAlpha Numeric value specifying the new alpha transparency level
#' (0-255) (Default: 80).
#'
#' @details
#' Not exported. This function takes a color or vector of colors in
#' hexadecimal or named colors and modifies their transparency by changing the
#' alpha channel value. It returns the modified color(s) with the updated
#' transparency.
#'
#' @return A character vector of modified colors with adjusted transparency in
#' hexadecimal.
#'
#' @examples
#' # Change transparency of a color
#' newColor <- changeColorAlpha(colorVec = "#FF0000", newAlpha = 80)
#'
#' @author Brian S. Roberts
#'
changeColorAlpha <- function(colorVec, newAlpha = 80) {
    transColX <- unlist(lapply(colorVec,function(charFillColX){
        rgbX <- as.numeric(col2rgb(charFillColX))
        transColX <- rgb(red= rgbX[1],green= rgbX[2],
            blue= rgbX[3],alpha= newAlpha, maxColorValue=255)
        return(transColX)
    }))
    return(transColX)
}

