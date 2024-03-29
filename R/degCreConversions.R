#' Calculate Raw Association Probability Odds Ratio (OR)
#'
#' Given a DegCre results list, this function calculates the raw association
#' probability odds ratio (OR) for each association.
#'
#' @param degCreResList List of DegCre results.
#'
#' @return A numeric vector of raw association probability odds ratios (OR)
#' for each association.
#'
#' @details
#' This function calculates the raw association probability odds ratio (OR)
#' for each association in a DegCre analysis output.
#' The OR is calculated relative to the distance bin null association
#' probability, which would happen if all CRE p-values were identical.
#' Thus it is a measure of the increase in association probability due to
#' CRE p-value information content over what would occur by random chance.
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load test data.
#' data(DexNR3C1)
#' 
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
#'                                    DegP=subDegGR$pVal,
#'                                    DegLfc=subDegGR$logFC,
#'                                    CreGR=subCreGR,
#'                                    CreP=subCreGR$pVal,
#'                                    CreLfc=subCreGR$logFC)
#'
#' #Calculate raw odds ratio.
#' ORvec <- calcRawAssocProbOR(degCreResListDexNR3C1)
#'
#' @author Brian S. Roberts
#'
#' @export
calcRawAssocProbOR <- function(degCreResList){

    hitsX <- degCreResList$degCreHits
    alphaValX <- degCreResList$alphaVal
    degPadjX <- S4Vectors::mcols(hitsX)$DegPadj
    binIdsX <- S4Vectors::mcols(hitsX)$distBinId
    rawAssocProbsX <- S4Vectors::mcols(hitsX)$rawAssocProb

    nullExpectedProbByBin <- tapply(degPadjX,binIdsX,function(degPadjY){
        expectTot <- (1-alphaValX)*length(which(degPadjY<=alphaValX))
        expectProbY <- expectTot/length(degPadjY)
        return(expectProbY)
    })
    uniqBinIds <- as.integer(names(nullExpectedProbByBin))
    nullExpectedProbByBin <- as.numeric(nullExpectedProbByBin)

    outRawAssocProbOR <- rep(NA,length(rawAssocProbsX))

    for(binIdY in uniqBinIds){
        maskY <- which(binIdsX == binIdY)
        maskBinY <- which(uniqBinIds == binIdY)

        nullExpectedProbBinY <- nullExpectedProbByBin[maskBinY]

        if(nullExpectedProbBinY == 0){
            outORY <- rep(1,length(maskY))
        }
        else{
            outORY <- rawAssocProbsX[maskY]/nullExpectedProbBinY
        }

        outRawAssocProbOR[maskY] <- outORY
    }
    return(outRawAssocProbOR)
}


#' Convert DegCre Results List to GInteractions Object
#'
#' Given a DegCre results list, this function converts it into a GInteractions
#' object.
#'
#' @param degCreResList List of DegCre results.
#' @param assocAlpha The significance threshold for associations to be included
#' in the output (Default: \code{0.05}).
#'
#' @return A \link[InteractionSet]{GInteractions} object containing the
#' associations that pass the specified significance threshold.
#'
#' The \link[InteractionSet]{GInteractions} object has same metadata columns as
#' the \link[S4Vectors]{Hits} returned from \link{runDegCre} with additional
#' columns.
#' These additional columns are every metadata column in the input \code{DegGR}
#' or \code{CreGR} proceeded by either \code{Deg_} or \code{Cre_} in the
#' colname.
#'
#' @details
#' This function takes a DegCre results list as input and extracts the
#' significant associations based on the  \code{assocProbFDR} compared to the
#'specified significance threshold \code{assocAlpha}. It then creates a
#' \link[InteractionSet]{GInteractions} object  metadata columns from the input
#' list.
#'
#' If no associations pass the significance threshold, the function returns
#' NA' and prints a message.
#'
#' @importFrom InteractionSet GInteractions
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load test data.
#' data(DexNR3C1)
#' 
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
#'                                    DegP=subDegGR$pVal,
#'                                    DegLfc=subDegGR$logFC,
#'                                    CreGR=subCreGR,
#'                                    CreP=subCreGR$pVal,
#'                                    CreLfc=subCreGR$logFC)
#'
#' #Create GInteractions object.
#' gInteractions <-
#'  convertdegCreResListToGInteraction(degCreResList=degCreResListDexNR3C1,
#'                                     assocAlpha = 0.01)
#'
#' @author Brian S. Roberts
#'
#' @export
convertdegCreResListToGInteraction <- function(degCreResList,
                                               assocAlpha = 0.05) {

    DegGRX <- degCreResList$DegGR
    CreGRX <- degCreResList$CreGR

    degCreHits <- degCreResList$degCreHits

    maskPassAlpha <- which(S4Vectors::mcols(degCreHits)$assocProbFDR<=
        assocAlpha)

    if(length(maskPassAlpha)<1){
        warning("No associations FDRs pass alpha.")
        outGInter <- NA
    }
    else{
        keepDegCreHits <- degCreHits[maskPassAlpha]

        outGInter <-
         InteractionSet::GInteractions(granges(DegGRX[S4Vectors::queryHits(keepDegCreHits)]),
                                       granges(CreGRX[S4Vectors::subjectHits(keepDegCreHits)]))

        mcolsHitsDf <- data.frame(S4Vectors::mcols(keepDegCreHits))

        mcolsDegDf <- 
          data.frame(S4Vectors::mcols(DegGRX[S4Vectors::queryHits(keepDegCreHits)]))
        colnames(mcolsDegDf) <- paste("Deg",colnames(mcolsDegDf),sep="_")

        mcolsCreDf <- 
          data.frame(S4Vectors::mcols(CreGRX[S4Vectors::subjectHits(keepDegCreHits)]))
        colnames(mcolsCreDf) <- paste("Cre",colnames(mcolsCreDf),sep="_")

        mcolsAllDf <- data.frame(mcolsHitsDf,mcolsDegDf,mcolsCreDf)
        S4Vectors::mcols(outGInter) <- mcolsAllDf
    }
    return(outGInter)
}


#' Convert DegCre Results List to DataFrame
#'
#' Given a DegCre results list, this function converts it into a DataFrame for
#' further analysis and export.
#'
#' @param degCreResList List of DegCre results.
#' @param assocAlpha The significance threshold for associations to be included
#' in the output (Default: \code{0.05}).
#'
#' @return A \link[S4Vectors]{DataFrame} containing the significant associations
#' that pass the specified significance threshold. It is roughly in BEDPE
#' format.
#'
#' @details
#' This function takes a DegCre results list as input and extracts the
#' significant associations based on the adjusted p-values \code{assocProbFDR}
#' compared to the specified significance threshold \code{assocAlpha}.
#' It then creates a \link[S4Vectors]{DataFrame} with the genomic coordinates
#' of the significant associations from both the \code{DegGR} and \code{CreGR}
#' components of the input list.
#' These are marked as \code{Deg_} or \code{Cre_} with \code{chr}, \code{start},
#' \code{end}, and \code{strand}.
#' The coordinates are followed by the metadata of the \link[S4Vectors]{Hits}
#' \link[S4Vectors]{DataFrame}d by \link{runDegCre}.
#' These are then followed by all metadata columns in the input \code{DegGR} or
#' \code{CreGR} proceeded by either \code{Deg_} or \code{Cre_} in the colname.
#'
#' If no associations pass the significance threshold,
#' the function returns \code{NA}.
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load test data.
#' data(DexNR3C1)
#' 
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
#'                                    DegP=subDegGR$pVal,
#'                                    DegLfc=subDegGR$logFC,
#'                                    CreGR=subCreGR,
#'                                    CreP=subCreGR$pVal,
#'                                    CreLfc=subCreGR$logFC)
#'
#' #Create DataFrame.
#' outDf <-
#'   convertDegCreDataFrame(degCreResList=degCreResListDexNR3C1,
#'                                     assocAlpha = 0.05)
#'
#' #Write out as text file.
#' degCreDfFile <- tempfile(pattern="myDegCreResults",fileext=".tsv")
#'
#' write.table(outDf,file=degCreDfFile[1],sep="\t",row.names=FALSE,quote=FALSE)
#'
#' unlink(degCreDfFile[1])
#'
#' @author Brian S. Roberts
#'
#' @export
convertDegCreDataFrame <- function(degCreResList,
                                   assocAlpha = 0.05){

    DegGRX <- degCreResList$DegGR
    CreGRX <- degCreResList$CreGR

    degCreHits <- degCreResList$degCreHits

    maskPassAlpha <- which(S4Vectors::mcols(degCreHits)$assocProbFDR<=
        assocAlpha)

    if(length(maskPassAlpha)<1){
        warning("No associations FDRs pass alpha.")
        outDf <- NA
    }
    else{
        keepDegCreHits <- degCreHits[maskPassAlpha]

        keepDegGRx <- DegGRX[S4Vectors::queryHits(keepDegCreHits)]
        keepDegDfx <- as.data.frame(keepDegGRx)
        #get rid of "width" column
        keepDegDfx <- keepDegDfx[,c(seq_len(3),
                        seq(from=5,to=ncol(keepDegDfx)))]
        colnames(keepDegDfx)[1] <- "chr"
        colnames(keepDegDfx) <- paste("Deg",colnames(keepDegDfx),sep="_")

        keepCreGRx <- CreGRX[S4Vectors::subjectHits(keepDegCreHits)]
        keepCreDfx <- as.data.frame(keepCreGRx)
        #get rid of "width" column
        keepCreDfx <- keepCreDfx[,c(seq_len(3),
                        seq(from=5,to=ncol(keepCreDfx)))]
        colnames(keepCreDfx)[1] <- "chr"
        colnames(keepCreDfx) <- paste("Cre",colnames(keepCreDfx),sep="_")

        hitsMcolsDf <- as.data.frame(S4Vectors::mcols(keepDegCreHits))

        #keep only those that ar not redundant with those in Deg or Cre Dfs
        keepHitsMcolsDf <- 
          hitsMcolsDf[,c(seq_len(4),seq(from=8,to=10)),drop=FALSE]

        outDf <- data.frame(keepDegDfx[,seq_len(4),drop=FALSE],
            keepCreDfx[,seq_len(4),drop=FALSE],
            keepHitsMcolsDf,
            keepDegDfx[,seq(from=5,to=ncol(keepDegDfx)),drop=FALSE],
            keepCreDfx[,seq(from=5,to=ncol(keepCreDfx)),drop=FALSE])
    }

    return(outDf)
}


