#' Get Expected Associations per DEG
#'
#' Calculates the expected associations per DEG (Differentially Expressed Gene).
#'
#' @param degCreResList A list of DegCre results.
#' @param geneNameColName Character value of the name of the column in DegGR
#' (that was inputted to \link{runDegCre}) that contains gene names. If NULL,
#' the function will attempt to automatically find the gene name column.
#' (Default: \code{NULL})
#' @param assocAlpha Numeric value from 0 to 1 specifying the significance
#' threshold for associations. (Default: \code{0.05})
#'
#' @details
#' This function calculates the expected associations per DEG based on DegCre
#' analysis results. It first filters significant associations based on the
#' provided association significance threshold (\code{assocAlpha}) and then
#' computes the expected associations per gene. The function returns a
#' \link[S4Vectors]{DataFrame} with gene-level information, including expected
#' associations, number of associations, and significance thresholds.
#'
#' @return A \link[S4Vectors]{DataFrame} with the all data in the input DegGR
#' with these columns added:
#' \describe{
#'   \item{geneName}{Character values of gene names extracted from
#'   \code{geneNameColName} column (or column found if
#'   \code{geneNameColName = NULL}) in DegGR.}
#'   \item{expectAssocs}{Numeric values of the expected associations per gene.}
#'   \item{nAssocs}{Integer value of the number of associations passing
#'   \code{assocAlpha}  per gene.}
#'   \item{assocAlpha}{Numeric value from 0 to 1 of input \code{assocAlpha}}
#'   \item{degAlpha}{Numeric value from 0 to 1 of the significance threshold
#'      for DEGs. Obtained from \code{degCreResList}}
#' }
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load example data.
#' data(DexNR3C1)
#'
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomicRanges::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomicRanges::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
#'                                    DegP=subDegGR$pVal,
#'                                    DegLfc=subDegGR$logFC,
#'                                    CreGR=subCreGR,
#'                                    CreP=subCreGR$pVal,
#'                                    CreLfc=subCreGR$logFC)
#'
#'
#' # Get expected associations per DEG
#' expectAssocsDf <- getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
#'                                        geneNameColName = "GeneSymb",
#'                                        assocAlpha = 0.05)
#'
#' head(expectAssocsDf)
#'
#'
#' @author Brian S. Roberts
#'
#' @export
getExpectAssocPerDEG <- function(degCreResList,
                                 geneNameColName = NULL,
                                 assocAlpha = 0.05) {

    degCreHits <- degCreResList$degCreHits
    DegGR <- degCreResList$DegGR
    alphaDeg <- degCreResList$alphaVal

    #get query level expected assocs
    maskPassAssocFDR <- which(mcols(degCreHits)$assocProbFDR <= assocAlpha)

    passDegCreHits <- degCreHits[maskPassAssocFDR]

    #now convert to gene level
    allColNameDegGr <- colnames(mcols(DegGR))

    if(is.null(geneNameColName)){
        maskGeneColName <- grep("gene",allColNameDegGr,ignore.case=TRUE)
        if(length(maskGeneColName)!=1){
            warning("Can not find unique gene column name in DegGR mcols.
                    Re-run specifying <geneNameColName>.")
            return(NA)
        }
    }
    else{
        maskGeneColName <- which(allColNameDegGr == geneNameColName)
    }

    allGeneNames <- mcols(DegGR)[,maskGeneColName]

    allPassAssocProbs <- mcols(passDegCreHits)$assocProb

    expectAssocPerGene <-
        tapply(allPassAssocProbs,
        INDEX=allGeneNames[queryHits(passDegCreHits)],
        FUN=function(expectX){
            return(sum(expectX))
        })

    nAssocPerGeneRle <- rle(sort(allGeneNames[queryHits(passDegCreHits)]))

    rawNAssocPerGene <- nAssocPerGeneRle$lengths

    maskToUniqAllGenes <-
        charmatch(names(expectAssocPerGene),nAssocPerGeneRle$values)

    nAssocPerGene <- rawNAssocPerGene[maskToUniqAllGenes]


    allDegGRMcolsDf <- mcols(DegGR)

    uniqAllGenes <- unique(allGeneNames)

    uniqGeneIndices <- unlist(lapply(uniqAllGenes,function(uniqGeneX){
        maskGeneX <- which(allDegGRMcolsDf[,maskGeneColName]==uniqGeneX)
        return(maskGeneX[1])
    }))

    uniqGeneDegMcolsDf <- allDegGRMcolsDf[uniqGeneIndices,]

    expectAssocsList <- lapply(uniqGeneDegMcolsDf[,maskGeneColName],
        function(geneNameX){

            maskGeneInExpect <- which(names(expectAssocPerGene) == geneNameX)
            if(length(maskGeneInExpect)==0){
                expectAssocX <- NA
                nAssocX <- NA
            }
            else{
                expectAssocX <- expectAssocPerGene[maskGeneInExpect]
                nAssocX <- nAssocPerGene[maskGeneInExpect]
            }
            return(c(expectAssocX,nAssocX))
        })

    expectAssocsAndNmat <- matrix(unlist(expectAssocsList),ncol=2,byrow=TRUE)
    expectAssocs <- expectAssocsAndNmat[,1]
    nAssocs <- expectAssocsAndNmat[,2]

    maskNaExpectedAssocs <- which(is.na(expectAssocs))
    subMaskNaIsDeg <-
        which(uniqGeneDegMcolsDf$pAdj[maskNaExpectedAssocs] <= alphaDeg)

    expectAssocs[maskNaExpectedAssocs[subMaskNaIsDeg]] <- 0

    nAssocs[maskNaExpectedAssocs[subMaskNaIsDeg]] <- 0

    uniqGeneDegMcolsDf$expectAssocs <- expectAssocs
    uniqGeneDegMcolsDf$nAssocs <- nAssocs

    uniqGeneDegMcolsDf <-
        uniqGeneDegMcolsDf[order((1/uniqGeneDegMcolsDf$expectAssocs),
            uniqGeneDegMcolsDf$pAdj,na.last=TRUE),]

    uniqGeneDegMcolsDf$assocAlpha <- assocAlpha
    uniqGeneDegMcolsDf$degAlpha <- alphaDeg

    return(uniqGeneDegMcolsDf)
}

#' Plot Histogram of Expected Associations per DEG
#'
#' Plots a histogram of the expected number of associations per DEG
#' (Differentially Expressed Gene) based on DegCre analysis.
#'
#' @param expectAssocPerDegDf \link[S4Vectors]{DataFrame} output of
#' \link{getExpectAssocPerDEG}.
#' @param barOutlineColor Color for the outline of the histogram bars.
#' (Default: \code{#88CCEE})
#' @param barFillColor Fill color for the histogram bars. If \code{NULL},
#' it will be derived from \code{barOutlineColor} with adjusted transparency.
#' @param extraText Logical, indicating whether additional text information
#' (Details) should be added to the plot.
#'
#' @details
#' This function generates a histogram of the expected number of associations
#' per DEG and optionally adds additional text information to the plot, such
#' as DEG FDR, association FDR, and the fraction of DEGs with at least one
#' association.
#' Plot displays a dashed line a value indicating the median expected
#' associations per DEG.
#'
#' @return Invisibly, the median expected associations per DEG.
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load example data.
#' data(DexNR3C1)
#'
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomicRanges::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomicRanges::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
#'                                    DegP=subDegGR$pVal,
#'                                    DegLfc=subDegGR$logFC,
#'                                    CreGR=subCreGR,
#'                                    CreP=subCreGR$pVal,
#'                                    CreLfc=subCreGR$logFC)
#'
#' # Generate data frame of expected associations per DEG
#' expectAssocPerDegDf <-
#'  getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
#'                       geneNameColName = "GeneSymb",
#'                       assocAlpha = 0.05)
#'
#' # Plot histogram of expected associations per DEG
#' medianExpAssocs <- plotExpectedAssocsPerDeg(expectAssocPerDegDf,
#'                                             barOutlineColor = "blue",
#'                                             extraText = TRUE)
#'
#' @author Brian S. Roberts
#'
#' @export
plotExpectedAssocsPerDeg <- function(expectAssocPerDegDf,
                                     barOutlineColor = "#88CCEE",
                                     barFillColor = NULL,
                                     extraText = FALSE) {

    if(is.null(barFillColor)){
        barFillColor <- changeColorAlpha(barOutlineColor,newAlpha=150)
    }

    rawExpectAssocs <- expectAssocPerDegDf$expectAssocs

    expectAssocPerDegDf <-
        expectAssocPerDegDf[which(!is.na(rawExpectAssocs)),]

    expectAssocs <- expectAssocPerDegDf$expectAssocs

    assocAlpha <- expectAssocPerDegDf$assocAlpha[1]
    degAlpha <- expectAssocPerDegDf$degAlpha[1]

    medianExpAssocs <- median(expectAssocs)

    #calculate the fraction of DEGs with at least one assoc
    totalDegs <- length(which(expectAssocPerDegDf$pAdj<=degAlpha))
    degsWithAtLeastOneAssoc <- length(which(expectAssocs>0))

    fractionDegsWithAssocs <- degsWithAtLeastOneAssoc/totalDegs

    histOuts <- hist(expectAssocs,breaks=20,col=barFillColor,
        border=barOutlineColor,xlab="Expect. Assocs. Per DEG",
        freq=TRUE,ylab="Num. DEGs",main="")

    plotTopY <- max(histOuts$counts)
    textCex <- par("cex.lab") - 0.1
    textOffset <- 0.2*textCex

    lines(x=c(medianExpAssocs,medianExpAssocs),y=c(0,max(histOuts$counts)),
        col="black",lwd=1.5,lty="dashed")

    par(xpd=NA)

    text(x=medianExpAssocs,y=plotTopY,
        label=round(medianExpAssocs,2),col="black",pos=3,
        cex=textCex,offset=textOffset)

    if(extraText){
        degFDRText <- bquote("DEG FDR" <= .(degAlpha))

        assocFDRText <- bquote("Assoc. FDR" <= .(assocAlpha))

        fracDegText <- paste("Frac. Deg =",
            round(fractionDegsWithAssocs,2))

        barMids <- histOuts$mids
        midRange <- max(barMids) - min(barMids)
        extraTextX <-  min(barMids) +0.6*midRange

        par(xpd=NA)

        text(x=extraTextX,y=plotTopY,label=degFDRText,offset=textOffset,
            pos=1,cex=textCex)

        text(x=extraTextX,y=plotTopY,label=assocFDRText,
            offset=textOffset + 1.1*textCex,pos=1,
            cex=textCex)

        text(x=extraTextX,y=plotTopY,label=fracDegText,
            offset=textOffset + 2.1*textCex,pos=1,
            cex=textCex)
    }
    invisible(medianExpAssocs)
}
