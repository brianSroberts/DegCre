#' DegCre
#'
#' Probabilistic association of DEGs to CREs from differential data.
#'
#' @docType package
#' @name DegCre
#' @aliases DegCre
#' @aliases DegCre-package
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand width
#' @importFrom IRanges findOverlaps distance start end reduce tile
#' @importFrom GenomicRanges granges GRanges GRangesList makeGRangesFromDataFrame
#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importFrom S4Vectors subjectHits queryHits mcols Hits
#' @importClassesFrom S4Vectors Hits
#' @importFrom stats median rnorm pbinom quantile p.adjust
#' @importFrom graphics plot lines points par axis layout polygon text hist
#' @importFrom grDevices colorRampPalette dev.size rgb col2rgb
#' @importFrom utils head tail
#' @importFrom plotgardener colorby plotPairsArches annoHeatmapLegend plotSignal annoYaxis plotGenomeLabel plotGenes plotText pageCreate
"_PACKAGE"

#' DegCre input data for examples.
#'
#' @name DexNR3C1
#' @docType data
#' @author Brian S. Roberts
#' @details This is a list with two slots: DegGR and CreGR.
#' This data was derived from work by McDowell et al. in which they generated
#' RNA-seq and ChIP-seq data by treating A549 cells with dexamethasone at
#' several time points. Specifically this is RNA-seq and NR3C1 ChIP-seq at four
#' hours versus control.
#' @format A named list with two slots: DegGR and CreGR.
#' \describe{
#' \item{DegGR}{\link[GenomicRanges]{GRanges} of RNA-seq data.
#' The coordinates reference TSS sites. It has the following mcols:
#' \describe{
#' \item{promGeneName}{\href{https://epd.expasy.org/epd}{EPDNew} promoter names}
#' \item{GeneSymb}{Gene symbols}
#' \item{GeneID}{Ensembl gene ids}
#' \item{baseMean}{baseMean values
#' from \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2}}
#' \item{logFC}{Log-2 fold-changes
#' from \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2}}
#' \item{pVal}{P-values
#' from \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2}}
#' \item{pAdj}{Adjusted p-values
#' from \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2}}
#' }
#' }
#' \item{CreGR}{\link[GenomicRanges]{GRanges} of differntial CRE data.
#' The coordinates reference signal regions. It has the following mcols:
#' \describe{
#' \item{logFC}{Log-2 fold-changes from
#' \href{https://bioconductor.org/packages/release/bioc/html/csaw.html}{csaw}}
#' \item{pVal}{P-values from
#' \href{https://bioconductor.org/packages/release/bioc/html/csaw.html}{csaw}}
#' \item{pAdj}{Adjusted p-values from
#' \href{https://bioconductor.org/packages/release/bioc/html/csaw.html}{csaw}}
#' }
#' }
#' }
#' @references \url{https://genome.cshlp.org/content/28/9/1272}
#' @keywords data
NULL
#'
#' Generate DegCre associations
#'
#' Create DEG to CRE associations from differential data.
#'
#' @param DegGR A \link[GenomicRanges]{GRanges} object of gene TSSs. Multiple
#' TSSs per gene are allowed.
#' @param DegP A numeric vector of differential expression p-values for genes
#' in \code{DegGR}.
#' @param DegLfc A numeric vector of log fold-change values of differential
#' expression for gene in \code{DegGR}. Required when
#' \code{reqEffectDirConcord = TRUE}. (Default: \code{NULL})
#' @param CreGR A \link[GenomicRanges]{GRanges} object of CRE regions.
#' @param CreP A numeric vector differential signal p-values for regions in
#' \code{CreGR}.
#' @param CreLfc A numeric vector log fold-change values of differential signal
#' for regions in \code{CreGR}. Required when
#' \code{reqEffectDirConcord = TRUE}. (Default: \code{NULL})
#' @param reqEffectDirConcord A logical whether to require concordance between
#' the effect direction between DEG and CRE differential values.
#' (Default: \code{TRUE})
#' @param padjMethod A character value indicating the method for p-value
#' adjustment. Do not change from default under most circumstances. Can be any
#' method name accepted by \link[stats]{p.adjust} (Default: \code{bonferroni})
#' @param maxDist An integer value specifying the maximum distance for
#' probability calculation of TSS to CRE associations. (Default: \code{1e6})
#' @param verbose A logical indicating whether to print messages of step
#' completion and algorithm results. (Default: \code{TRUE})
#' @param smallestTestBinSize An integer value specifying the size
#' (number of elements) of the smallest distance bin to be considered in the
#' optimization algorithm. (Default: \code{100})
#' @param fracMinKsMedianThresh A numeric value between 0 and 1 specifying the
#' optimization criterion for the distance bin size algorithm (See Details).
#' (Default: \code{0.2})
#' @param alphaVal A numeric value between 0 and 1 specifying the alpha value
#' for DEG significance. (Default: \code{0.01})
#' @param binNOverride An integer value specifying the number of elements per
#' distance bin. When specified, overrides distance bin size optimization
#' (Not recommended). (Default: \code{NULL})
#'
#' @details
#' The DegCre algorithm considers experimental data from a perturbation
#' experiment and produces associations between cis-regulatory elements
#' (CREs) and differentially expressed genes (DEGs).
#' The user provides differential expression data such as RNA-seq, and
#' differential regulatory signal data such as ATAC-seq, DNase
#' Hypersensitivity, and ChIP-seq.
#' For RNA-seq analysis, we suggest methods such as
#' \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2}
#' or \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR}.
#' For the analysis of differential regulatory data we recommend
#' \href{https://bioconductor.org/packages/release/bioc/html/csaw.html}{csaw}.
#' As an example experiment, we use data from McDowell et al. (PMID = 30097539)
#' in which A549 cells were treated with dexamethasone and control.
#' RNA-seq and ChIP-seq data were collected at various time points.
#'
#' A complete description of the mathematical basis of the DegCre core
#' algorithms is provided in
#' \href{https://www.biorxiv.org/content/10.1101/2023.10.04.560923v1}{DegCre bioRxiv}.
#' DegCre takes two inputs. The first is a GRanges of p-values and optionally
#' log fold-changes associated with DEG TSSs.
#' The second input is a GRanges of differential signal p-values and optionally
#' log fold-changes for CRE regions.
#' DegCre generates a \link[S4Vectors]{Hits} object of all associations between
#' DEG TSSs and CREs within \code{maxDist}.
#' Associations are then binned by TSS-to-CRE distance according to an
#' algorithm that balances resolution (many bins with few members)
#' versus minimization of the deviance of each bin's CRE p-value distribution
#' from the global distribution, seleting an optimal bin size.
#'
#' Next, DegCre applies a non-parametric algorithm to find concordance between
#' and CRE differential effects within bins and derives an association
#' probability.
#' For all association probabilities involving one given CRE, the probabilities
#' are adjusted to favor associations across shorter distances.
#' An FDR of the association probability is then estimated. Results are
#' returned in list containing a \link[S4Vectors]{Hits} object and both
#' input GRanges.
#'
#' @return A named list containing:
#' \describe{
#'   \item{degCreHits}{A \link[S4Vectors]{Hits} object with metadata.
#'   The \link[S4Vectors]{queryHits} of
#'   \code{degCreHits} reference \code{DegGR}.
#'   The \link[S4Vectors]{subjectHits} of
#'   \code{degCreHits} reference \code{CreGR}}
#'   \item{binHeurOutputs}{List of outputs from the distance binning algorithm.}
#'   \item{alphaVal}{Numeric alpha value used for DEG significance threshold.}
#'   \item{DegGR}{\link[GenomicRanges]{GRanges} of input \code{DegGR} with
#'   added metadata columns "pVal", "pAdj",and possibly "logFC"
#'   if \code{reqEffectDirConcord==TRUE}. Will overwrite existing metadata
#'   with same colnames.}
#'   \item{CreGR}{\link[GenomicRanges]{GRanges} of input \code{CreGR} with
#'   added metadata columns "pVal", "pAdj",and possibly "logFC" if
#'   \code{reqEffectDirConcord==TRUE}. Will overwrite existing metadata with
#'   same colnames.}
#' }
#' The degCreHits \link[S4Vectors]{Hits} object metadata has these columns:
#' \describe{
#'     \item{assocDist}{Integer of distance in base pairs between the TSS and
#'     CRE for the association.}
#'     \item{assocProb}{Numeric from 0 to 1 of association probability.}
#'     \item{assocProbFDR}{Numeric from 0 to 1 of False discovery rate of
#'     the association probability exceeding distance only null.}
#'     \item{rawAssocProb}{Numeric from 0 to 1 of association probability not
#'     adjusted for DEG significance or shorter associations involving
#'     this CRE.}
#'     \item{CreP}{Numeric of differential p-value of the CRE.}
#'     \item{DegP}{Numeric of differential p-value of the DEG.}
#'     \item{DegPadj}{Numeric of differential adjusted p-value of the DEG.}
#'     \item{binAssocDist}{Integer of the maximum association distance cutoff
#'     for the bin containing the association.}
#'     \item{numObs}{Integer number of associations in the distance bin
#'     containing the association.}
#'     \item{distBinId}{Integer that uniquely identifies the distance
#'     containing the association.}
#' }
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load sample data.
#' data(DexNR3C1)
#' 
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' #With defaults.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
#'                                    DegP=subDegGR$pVal,
#'                                    DegLfc=subDegGR$logFC,
#'                                    CreGR=subCreGR,
#'                                    CreP=subCreGR$pVal,
#'                                    CreLfc=subCreGR$logFC)
#'
#' #With custom settings.
#' modDegCreResList <- runDegCre(DegGR=subDegGR,
#'                            DegP=subDegGR$pVal,
#'                            CreGR=subCreGR,
#'                            CreP=subCreGR$pVal,
#'                            reqEffectDirConcord=FALSE,
#'                            maxDist=1e5,
#'                            alphaVal=0.001)
#'
#' @author Brian S. Roberts
#'
#' @export
runDegCre <- function(DegGR,
                   DegP,
                   DegLfc = NULL,
                   CreGR,
                   CreP,
                   CreLfc = NULL,
                   reqEffectDirConcord = TRUE,
                   padjMethod = "bonferroni",
                   maxDist = 1e6,
                   verbose = TRUE,
                   smallestTestBinSize = 100,
                   fracMinKsMedianThresh = 0.2,
                   alphaVal = 0.01,
                   binNOverride = NULL){

  #check if padjMethod is TRUE and if log fold changes have been provided
    if(all(c(reqEffectDirConcord,any(c(is.null(DegLfc),is.null(CreLfc)))))){
        warning("Log fold changes not given when
              reqEffectDirConcord = TRUE; returning NA")
        return(NA)
    }

    outPadj <- rep(NA,length(DegP))
    maskNoNAPProm <- which(!is.na(DegP))
    noNaPadj <- p.adjust(DegP[maskNoNAPProm],method=padjMethod)
    outPadj[maskNoNAPProm] <- noNaPadj
    DegPadj <- outPadj

    #add pVals and pAdj to input GRs with controlled column names
    #for downstream functions
    DegGR$pVal <- DegP
    DegGR$pAdj <- DegPadj
    if(!is.null(DegLfc)){
        DegGR$logFC <- DegLfc
    }

    CreGR$pVal <- CreP
    if(!is.null(CreLfc)){
        CreGR$logFC <- CreLfc
    }

    #compare only those regions with non NA p-values in
    # both DegP, CreP, CrePadj, and DegPadj

    maskNoNaCreP <- which(!is.na(CreP))

    maskNoNaPromP <- which(!is.na(DegP))
    maskNoNaPromPadj <- which(!is.na(DegPadj))
    maskNoNaPromBoth <- intersect(maskNoNaPromP,maskNoNaPromPadj)

    origCreGR <- CreGR
    origDegGR <- DegGR

    origDegPadj <- DegPadj
    origCreP <- CreP

    CreGR <- CreGR[maskNoNaCreP]
    DegGR <- DegGR[maskNoNaPromBoth]

    CreP <- CreP[maskNoNaCreP]

    DegP <- DegP[maskNoNaPromBoth]
    DegPadj <- DegPadj[maskNoNaPromBoth]

    #get hits with distance via getAssocDistHits
    degCreHits <- getAssocDistHits(DegGR=DegGR,
        CreGR=CreGR,
        maxDist=maxDist)

    #if reqEffectDirConcord = TRUE, only those hits with log fold changes
    #in the same directions are considered in the remaining analysis
    if(reqEffectDirConcord){
        DegLfc <- DegLfc[maskNoNaPromBoth]
        CreLfc <- CreLfc[maskNoNaCreP]

        promHitSpaceLfcs <- DegLfc[S4Vectors::queryHits(degCreHits)]
        creHitSpaceLfcs <- CreLfc[S4Vectors::subjectHits(degCreHits)]
        maskConcordantHits <- which(sign(promHitSpaceLfcs)==
            sign(creHitSpaceLfcs))
        degCreHits <- degCreHits[maskConcordantHits]
    }

    #Want to bin hits by distance (from DegGR to CreGR).
    #Want bins to have roughly the same distribution of CreP values.
    #distBinHeuristic finds a bin size that balances
    #the deviation in CreP distribution from the whole against
    #having a large number of bins so that the resolution is high.
    #All bins will have the same number of observations except the last bin
    #which will almost always have more.

    if(is.null(binNOverride)){
        binHeuristicList <- distBinHeuristic(degCreHits=degCreHits,
        CreP=CreP,
        smallestTestBinSize=smallestTestBinSize,
        fracMinKsMedianThresh=fracMinKsMedianThresh,
        verbose=verbose)

        pickedBinSize <- binHeuristicList$pickedBinSize
    }
    else{
        binHeuristicList <- list()
        binHeuristicList$pickedBinSize <- binNOverride
        pickedBinSize <- binNOverride
    }

    if(verbose){
        message("Hit distance bin n = ",pickedBinSize)
    }
    #now bin data by distance into bins with pickedBinSize hits per bin

    #going to sort degCreHits by distance rather than queryHits which
    #makes Hits mad. Convert to a dataframe for a while
    promCreHitsDf <- as.data.frame(degCreHits)

    sortPromCreHitsDf <- 
      promCreHitsDf[order(promCreHitsDf$assocDist),,drop=FALSE]

    numHits <- nrow(sortPromCreHitsDf)

    allSortHitIndices <- seq_len(numHits)

    rawSortHitChunksList <- split(allSortHitIndices,
        ceiling(seq_along(allSortHitIndices)/(pickedBinSize)))

    #last bin can be small, fix if too small.
    #Too small is less than half pickedBinSize
    if(length(rawSortHitChunksList[length(rawSortHitChunksList)]) <
        0.8*pickedBinSize){
        #if last bin is too small, split the last two bins into
        #bins with equal number of hits
        last2RawSortChunkIndices <-
            c(rawSortHitChunksList[[(length(rawSortHitChunksList)-1)]],
            rawSortHitChunksList[[length(rawSortHitChunksList)]])

        last2MidPoint <- floor(length(last2RawSortChunkIndices)/2)

        rawSortHitChunksList[[(length(rawSortHitChunksList)-1)]] <-
            last2RawSortChunkIndices[seq_len(last2MidPoint)]

        rawSortHitChunksList[[length(rawSortHitChunksList)]] <-
            last2RawSortChunkIndices[(last2MidPoint+1):
                length(last2RawSortChunkIndices)]

        sortHitChunksList <- rawSortHitChunksList

        if(verbose){
            message("Small last bin resolved.")
        }
    }
    else{
        sortHitChunksList <- rawSortHitChunksList
    }


    #calculate association probs
    if(verbose){
        message("Running over ",
            length(sortHitChunksList)," hit distance thresholds")
    }

    #add bin membership data to HitsDf

    distBinId <- unlist(lapply(seq_along(sortHitChunksList),function(binIdX){
        numMembers <- length(sortHitChunksList[[binIdX]])
        return(rep(binIdX,numMembers))
    }))

    sortPromCreHitsDf <- data.frame(sortPromCreHitsDf,distBinId)

    #distance binning is complete. Now run probability calculations

    listStatsPerChunk <- lapply(sortHitChunksList,function(chunkIx){

        #first calc the probs of the CREs associating with a DEG
        chunkCREStatsMatX <- calcDependIndependEnrichStats(
            hitsWithDistDf=sortPromCreHitsDf,
            subHitsIndex=chunkIx,
            dependPadj=DegPadj,
            independP=CreP,
            alpha=alphaVal)

        maxAssocDist <- max(sortPromCreHitsDf$assocDist[chunkIx])

        outStatsMatX <- cbind(chunkCREStatsMatX,
            rep(maxAssocDist,nrow(chunkCREStatsMatX)))

        return(outStatsMatX)
    })

    rawAllDistBinsStatsMat <- do.call(rbind,listStatsPerChunk)
    #colnames of rawAllDistBinsStatsMat:
    #"independP","assocProb","totalObs","binAssocDist"

    #adjust association probs by the whether the DEG is
    #differentially expressed, i.e. if DEG padj <= alphaVal
    sortDegPadj <- DegPadj[sortPromCreHitsDf[,1]]

    sortDegTruePosProb <- rep(0,length(sortDegPadj))
    sortDegTruePosProb[which(sortDegPadj<=alphaVal)] <- 1

    adjAssocProb <- sortDegTruePosProb*rawAllDistBinsStatsMat[,2]

    sortDegP <- DegP[sortPromCreHitsDf[,1]]

    if(verbose){
        message("Adjusting association probability for distance.")
    }

    #correct adjAssocProb for distance
    corrAdjAssocProb <- correctAssocProbs(sortHitsDf=sortPromCreHitsDf,
        assocProbs=adjAssocProb,
        refAssocProbs=adjAssocProb)

    #collect results up to here
    allDistBinsStatsMat <- cbind(corrAdjAssocProb,rawAllDistBinsStatsMat[,2],
        rawAllDistBinsStatsMat[,1],sortDegP,sortDegPadj,
        rawAllDistBinsStatsMat[,4],rawAllDistBinsStatsMat[,3])

    colnames(allDistBinsStatsMat) <- c("assocProb","rawAssocProb",
    "CreP","DegP","DegPadj","binAssocDist","numObs")


    if(verbose){
        message("Calculation of assoc probabilities over hit bins complete.")
    }

    #calculate FDR from binomial distribution
    listBinomStatsPerChunk <- lapply(sortHitChunksList,
         function(chunkIx){

         assocProbFDRX <-
         calcBinomFDRperBin(allDistBinsStatsMat=allDistBinsStatsMat,
            chunkI=chunkIx,
            alphaVal=alphaVal)

         return(assocProbFDRX)
    })

    allAssocProbFDR <- unlist(listBinomStatsPerChunk)


    if(verbose){
        message("Assoc FDR calculation complete.")
    }


    #add Ps  of associations; reorder too
    allDistBinsStatsMat <- cbind(allDistBinsStatsMat[,1,drop=FALSE],
        allAssocProbFDR,
        allDistBinsStatsMat[,2:ncol(allDistBinsStatsMat),drop=FALSE])

    colnames(allDistBinsStatsMat)[1] <- "assocProb"
    colnames(allDistBinsStatsMat)[2] <- "assocProbFDR"

    #put together allDistBinsStatsMat and sortPromCreHitsDf
    sortPromCreHitsStatsDf <- data.frame(sortPromCreHitsDf[,seq_len(3),
                                                           drop=FALSE],
        allDistBinsStatsMat,sortPromCreHitsDf[,4,drop=FALSE])
    
    colnames(sortPromCreHitsStatsDf)[ncol(sortPromCreHitsStatsDf)] <-
        "distBinId"

    #query sort sortPromCreHitsStatsDf
    sortPromCreHitsStatsDf <-
        sortPromCreHitsStatsDf[order(sortPromCreHitsStatsDf[,1]),]

    #now convert the indices of sortPromCreHitsStatsDf back into the
    #indices of origCreGR and origDegGR
    origDegGRIndices <- maskNoNaPromBoth[sortPromCreHitsStatsDf[,1]]
    origCreGRIndices <- maskNoNaCreP[sortPromCreHitsStatsDf[,2]]

    sortPromCreHitsStatsDf[,1] <- origDegGRIndices
    sortPromCreHitsStatsDf[,2] <- origCreGRIndices

    #convert sortPromCreHitsStatsDf back into a Hits object
    sortPromCreHitsStats <- S4Vectors::Hits(from=sortPromCreHitsStatsDf[,1],
        to=sortPromCreHitsStatsDf[,2],
        nLnode=length(origDegGR),
        nRnode=length(origCreGR),
        sort.by.query=TRUE)


    mcolsHitsDf <- 
      sortPromCreHitsStatsDf[,3:ncol(sortPromCreHitsStatsDf),drop=FALSE]

    #add metadata columns to sortPromCreHitsStatsDf
    S4Vectors::mcols(sortPromCreHitsStats) <- mcolsHitsDf

    #form output list
    outList <- list()
    outList$degCreHits <- sortPromCreHitsStats
    outList$binHeurOutputs <- binHeuristicList
    outList$alphaVal <- alphaVal
    outList$DegGR <- origDegGR
    outList$CreGR <- origCreGR

    if(verbose){
        message("DegCre calculations complete.")
    }

    return(outList)
}

#' Run DegCre with DEG alpha optimization.
#'
#' Runs DegCre across a set of DEG alpha thesholds to find optimal performance.
#'
#' @param DegGR A \link[GenomicRanges]{GRanges} object of gene TSSs. Multiple
#' TSSs per gene are allowed.
#' @param DegP A numeric vector of differential expression p-values for genes
#' in \code{DegGR}.
#' @param DegLfc A numeric vector of log fold-change values of differential
#' expression for gene in \code{DegGR}. Required when
#' \code{reqEffectDirConcord = TRUE}. (Default: \code{NULL})
#' @param CreGR A \link[GenomicRanges]{GRanges} object of CRE regions.
#' @param CreP A numeric vector differential signal p-values for regions in
#' \code{CreGR}.
#' @param CreLfc A numeric vector log fold-change values of differential
#' signal for regions in \code{CreGR}. Required when
#' \code{reqEffectDirConcord = TRUE}. (Default: \code{NULL})
#' @param reqEffectDirConcord A logical whether to require concordance between
#' the effect direction between DEG and CRE differential values.
#' (Default: \code{NULL})
#' @param padjMethod A character value indicating the method for p-value
#' adjustment. Do not change from default under most circumstances. Can be any
#' method name accepted by \code{p.adjust()} (Default: \code{bonferroni})
#' @param maxDist An integer value specifying the maximum distance for
#' probability calculation of TSS to CRE associations. (Default: \code{1e6})
#' @param verbose A logical indicating whether to print messages of step
#' completion and algorithm results. (Default: \code{NULL})
#' @param smallestTestBinSize An integer value specifying the size
#' (number of elements) of the smallest distance bin to be considered in the
#' optimization algorithm. (Default: \code{100})
#' @param fracMinKsMedianThresh A numeric value between 0 and 1 specifying the
#' optimization criterion for the distance bin size algorithm (See Details).
#' (Default: \code{0.2})
#' @param testedAlphaVals A numeric vector of DEG alpha values to test
#' (Default: \code{c(0.005,0.01,0.02,0.03,0.05,0.1)}).
#' @param minNDegs An integer specifying minimum number of DEGs that pass
#' the lowest \code{testedAlphaVals}. (Default: \code{5})
#'
#' @return A named list containing:
#' \describe{
#'   \item{alphaPRMat}{A matrix of Precision-Recall Area Under the Curve
#'   (AUC) values.}
#'   \item{degCreResListsByAlpha}{Named list of DegCre results lists indexed
#'   by the \code{testedAlphaVals}}.
#' }
#'
#' The columns of \code{alphaPRMat} are:
#' \describe{
#'   \item{alphaVal}{Numeric vector of tested DEG alpha value.}
#'   \item{AUC}{Numeric vector of Area under the curve of a Precision-Recall
#'   (PR) curve based on associations recovering significant DEGs.}
#'   \item{deltaAUC}{Numeric vector of PR AUC minus the AUC of the
#'   no-skill line.}
#'   \item{normDeltaAUC}{Numeric vector of the value of \code{deltaAUC}
#'   divided by one minus the no-skill AUC.}
#' }
#'
#' @details
#' This function runs \link{runDegCre} for each value in \code{testedAlphaVals}.
#' The performance at each tested alpha is evaluated with \link{degCrePRAUC}.
#' which generates a Precision-Recall curve based on the recovery rate of DEGs
#' by associations.
#' Various AUCs are calculated as performance metrics. Using the alpha with
#' the highest value of
#' \code{normDeltaAUC} is recommended (see Examples).

#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load sample data.
#' data(DexNR3C1)
#' 
#' subDegGR <-
#'  DexNR3C1$DegGR[which(GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1")]
#' subCreGR <-
#'  DexNR3C1$CreGR[which(GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1")]
#'
#' # Run DegCre over range of alpha values:
#' alphaOptList <- optimizeAlphaDegCre(DegGR = subDegGR,
#'                            DegP = subDegGR$pVal,
#'                            DegLfc = subDegGR$logFC,
#'                            CreGR = subCreGR,
#'                            CreP = subCreGR$pVal,
#'                            CreLfc = subCreGR$logFC)
#'
#' bestAlphaId <- which.max(alphaOptList$alphaPRMat[,4])
#' bestDegCreResList <- alphaOptList$degCreResListsByAlpha[[bestAlphaId]]
#'
#' @author Brian S. Roberts
#'
#' @export
optimizeAlphaDegCre <- function(DegGR,
                                DegP,
                                DegLfc=NULL,
                                CreGR,
                                CreP,
                                CreLfc=NULL,
                                reqEffectDirConcord=TRUE,
                                padjMethod="bonferroni",
                                maxDist=1e6,
                                verbose=FALSE,
                                smallestTestBinSize=100,
                                fracMinKsMedianThresh=0.2,
                                testedAlphaVals=c(0.005,0.01,0.02,
                                0.03,0.05,0.1),
                                minNDegs=5){
    
    #do not test DEG alpha vals so low that it does not result in at least
    #minNDegs passing alpha. The optimization algorithm requires a quorum
    #of passing DEGs to be valid
    tempPadjs <- p.adjust(DegP,method = padjMethod)
    
    nDegsPass <- unlist(lapply(testedAlphaVals,function(alphaX){
      return(length(which(tempPadjs<=alphaX)))
    }))
    
    maskAlphaPass <- which(nDegsPass>=minNDegs)
    
    if(length(maskAlphaPass)<length(testedAlphaVals) & length(maskAlphaPass)>0){
      testedAlphaVals <- testedAlphaVals[maskAlphaPass]
      testAlphaWarn <- 
        paste("Dropping tested alphas resulting in too few DEGS.",
              paste("New tested alphas =",paste(testedAlphaVals,collapse=",")),
              sep="\n")
      warning(testAlphaWarn)
    }
    
    #check if none of the testedAlphaVals meet too few DEG criteria
    #if do not run optimization and force a pickedAlpha
    runOptim <- TRUE
    
    if(length(maskAlphaPass)==0){
      #how many values are non 1?
      maskNon1 <- which(tempPadjs<1)
      nNon1 <- length(maskNon1)
      
      if(length(nNon1)<minNDegs){
        #there are not enough non-zero pAdjs to test.
        #Bypass opimization and set optimal alpha to max non-zero
        runOptim <- FALSE
        pickedAlpha <- max(tempPadjs[maskNon1])
        outMat <- matrix(c(pickedAlpha,0,0,0),nrow=1)
        
        outList <- list()
        outList$alphaPRMat <- outMat
        outList$degCreResListsByAlpha <- NA
        
        warnMsg <- 
          "Not enough non unity DEG adjusted pvalues. Returning highest non-unity"
        
        warning(warnMsg)
      }
      else{
        non1Padjs <- tempPadjs[maskNon1]
        minQProb <- minNDegs/length(non1Padjs)
        maxQProb <- nNon1/length(non1Padjs)
        midQProb <- mean(c(minQProb,maxQProb))
        
        testProbs <- c(minQProb,midQProb,maxQProb)
        
        testedAlphaVals <- quantile(tempPadjs,probs=testProbs)
        testAlphaWarn <- 
          paste("Recalculated alphas to make usable DEG sets.",
                paste("New tested alphas =",
                      paste(testedAlphaVals,collapse=",")),
                sep="\n")
        warning(testAlphaWarn)
      }
    }
    
    if(runOptim){
      
      alphaValNames <- paste("alpha",testedAlphaVals,sep="_")
      names(testedAlphaVals) <- alphaValNames
      
      degResForDistBinN <- runDegCre(DegGR=DegGR,
                                     DegP=DegP,
                                     DegLfc=DegLfc,
                                     CreGR=CreGR,
                                     CreP=CreP,
                                     CreLfc=CreLfc,
                                     reqEffectDirConcord=reqEffectDirConcord,
                                     padjMethod="bonferroni",
                                     maxDist=maxDist,
                                     verbose=verbose,
                                     alphaVal=0.01,
                                     fracMinKsMedianThresh=fracMinKsMedianThresh,
                                     smallestTestBinSize=smallestTestBinSize)
      
      pickedDistBinN <- degResForDistBinN$binHeurOutputs$pickedBinSize
      
      listByAlpha <- lapply(testedAlphaVals,function(alphaValX){
        if(verbose){
          message("Testing alpha = ",alphaValX)
        }
        
        degCreOutX <- runDegCre(DegGR=DegGR,
                                DegP=DegP,
                                DegLfc=DegLfc,
                                CreGR=CreGR,
                                CreP=CreP,
                                CreLfc=CreLfc,
                                reqEffectDirConcord=reqEffectDirConcord,
                                padjMethod="bonferroni",
                                maxDist=maxDist,
                                verbose=verbose,
                                alphaVal=alphaValX,
                                fracMinKsMedianThresh=fracMinKsMedianThresh,
                                binNOverride=pickedDistBinN,
                                smallestTestBinSize=smallestTestBinSize)
        
        #now get PR Auc values
        #6/27/2024
        #try a different optimization strategy
        #this approach is to maxmize the curve of unique CRE pVals
        #vs RawAssocProbOR
        
        #get rawORs
        rawORs <- calcRawAssocProbOR(degCreOutX)
        
        rawORsPerUniqCrePval <- 
          tapply(rawORs,
                 INDEX=mcols(degCreOutX$degCreHits)$CreP,
                 FUN=median)
        
        rawORsPerUniqCrePval <- as.numeric(rawORsPerUniqCrePval)
        uniqCrePs <- as.numeric(names(rawORsPerUniqCrePval))
        
        #order by pvals
        rawORsPerUniqCrePval <- 
          rawORsPerUniqCrePval[order(uniqCrePs,decreasing=TRUE)]
        
        #get AUC
        rawORAUC <- calcAUC(xVals=c(1:length(rawORsPerUniqCrePval)),
                            yVals=rawORsPerUniqCrePval)
        
        # PRAucResX <- degCrePRAUC(degCreOutX,
        #                          makePlot=FALSE,
        #                          nShuff=10,
        #                          alphaVal=alphaValX)
        
        outList <- degCreOutX
        outList$binHeurOutputs <- degResForDistBinN$binHeurOutputs
        outList$AUC <- rawORAUC
        #outList$deltaAUC <- PRAucResX$deltaAUC
        #outList$normDeltaAUC <- PRAucResX$normDeltaAUC
        return(outList)
      })
      
      allAUC <- unlist(lapply(listByAlpha,function(x){
        return(x$AUC)
      }))
      
      # allDeltaAUC <- unlist(lapply(listByAlpha,function(x){
      #   return(x$deltaAUC)
      # }))
      # 
      # allNormDeltaAUC <- unlist(lapply(listByAlpha,function(x){
      #   return(x$normDeltaAUC)
      # }))
      
      outMat <- cbind(testedAlphaVals,allAUC)
      colnames(outMat) <- c("alphaVal","AUC")
      
      # outMat <- cbind(testedAlphaVals,allAUC,allDeltaAUC,allNormDeltaAUC)
      # colnames(outMat) <- c("alphaVal","AUC","deltaAUC","normDeltaAUC")
      
      outList <- list()
      
      #clean up listByAlpha to not have AUC values
      
      listKeepByAlpha <- lapply(listByAlpha,function(listX){
        outList <- listX[seq_len(5)]
        return(outList)
      })
      
      outList$alphaPRMat <- outMat
      outList$degCreResListsByAlpha <- listKeepByAlpha
    }
    
    return(outList)
}

#' Calculate PR AUC for DegCre results.
#'
#' This function calculates the Precision-Recall Area Under the Curve (AUC)
#' from a DegCre results list.
#'
#' @param degCreResList A list of DegCre results.
#' @param makePlot Logical indicating whether to generate a plot of the
#' Precision-Recall curve. (Default: \code{TRUE})
#' @param nShuff Integer number of shuffles for no-skill curve.
#' (Default: \code{100})
#' @param alphaVal Numeric from 0 to 1 threshold alpha value of DEG
#' significance. (Default: alpha value from \code{degCreResList})
#' @param nThresh Integer number of threshold values for the Precision-Recall
#' curve. (Default: \code{200})
#'
#' @return Invisibly, a list containing:
#' \describe{
#'   \item{actualTprPpvMat}{A matrix of actual True Positive Rate (TPR) and
#'   apparent Positive Predictive Value (PPV).}
#'   \item{shuffTprQMat}{A matrix of shuffled TPR quantiles.}
#'   \item{shuffPpvQMat}{A matrix of shuffled PPV quantiles.}
#'   \item{AUC}{Numeric of the total Area Under the Curve (AUC) for the
#'   Precision-Recall curve.}
#'   \item{deltaAUC}{Numeric of the difference in AUC between the actual curve
#'   and shuffled curves.}
#'   \item{normDeltaAUC}{Numeric of the normalized difference in AUC.}
#' }
#'
#' @details
#' This function calculates the Precision-Recall curve and AUC based on the
#' provided DegCre results. It also estimates the statistical significance of
#' the AUC by shuffling the associations and calculating AUC for shuffled data.
#' Note that the PR AUCs tend to be small (0.05-0.2). Under the calculation
#' framework, a PR AUC of 1 could only be achieved from DegCre results in
#' which every association involves a significant DEG and has an association
#' probability of 1.
#' This situation will never actually occur but serves as a theoretical
#' optimum for comparison.
#'
#'
#' @examples
#' #Load required packages.
#' library(GenomicRanges)
#' 
#' #Load sample data.
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
#' #Plot PR curve.
#'
#' degCrePRAUC(degCreResList=degCreResListDexNR3C1)
#'
#' #Get PR results with out plotting.
#'
#' prAUCList <- degCrePRAUC(degCreResList=degCreResListDexNR3C1,
#'                         makePlot=FALSE)
#'
#' @author Brian S. Roberts
#'
#' @export
degCrePRAUC <- function(degCreResList,
                        makePlot=TRUE,
                        nShuff=100,
                        alphaVal=degCreResList$alphaVal,
                        nThresh=200){

       hitsDegCre <- degCreResList$degCreHits
    assocProb <- S4Vectors::mcols(hitsDegCre)$assocProb

    expectedDEGPos <- rep(0,length(hitsDegCre))
    expectedDEGPos[which(S4Vectors::mcols(hitsDegCre)$DegPadj<=alphaVal)] <- (1-alphaVal)
    totExpectDEGPos <- sum(expectedDEGPos)

    #vary the non-zero assocProb
    maskAdjProbNot0 <- which(assocProb>0)
    testedAdjProb <- assocProb[maskAdjProbNot0]
    testedExpectedDEGPos <- expectedDEGPos[maskAdjProbNot0]

    sortTestedAdjProb <- sort(unique(testedAdjProb),decreasing=TRUE)

    if(length(sortTestedAdjProb)<=nThresh){
        threshAdjProbs <- sortTestedAdjProb
    }
    else{
        threshAdjProbSet1 <- sortTestedAdjProb[seq_len(floor(nThresh/2))]
        maskSet2 <- c((floor(nThresh/2)+1):length(sortTestedAdjProb))
        qProbs <- rev(seq(from=0,to=1,length.out=nThresh-floor(nThresh/2)))
        threshAdjProbSet2 <- quantile(sortTestedAdjProb[maskSet2],qProbs)
        threshAdjProbs <- c(threshAdjProbSet1,threshAdjProbSet2)
    }

    listActualTPrPPV <- lapply(threshAdjProbs,
        function(threshAdjProbX){

        maskPassX <- which(testedAdjProb>=threshAdjProbX)
        expectDegAssocs <- sum(testedAdjProb[maskPassX]*
            testedExpectedDEGPos[maskPassX])
        tprX <- expectDegAssocs/totExpectDEGPos

        ppvX <- sum(testedAdjProb[maskPassX])/length(maskPassX)
        return(c(tprX,ppvX))
    })

    actualTPrPPVMAt <- matrix(unlist(listActualTPrPPV),ncol=2,byrow=TRUE)
    colnames(actualTPrPPVMAt) <- c("TPR","PPV")

    listShuffTprPPV <- lapply(seq_len(nShuff),
        function(nShuffi){

        shuffAdjAssocProbs <- sample(testedAdjProb,length(testedAdjProb),
            replace=FALSE)

        listShuffTprPpvI <- unlist(lapply(threshAdjProbs,
            function(threshAdjProbX){

            maskPassX <- which(testedAdjProb>=threshAdjProbX)
            expectDegAssocs <- sum(shuffAdjAssocProbs[maskPassX]*
                testedExpectedDEGPos[maskPassX])
            tprX <- expectDegAssocs/totExpectDEGPos

            ppvX <- sum(shuffAdjAssocProbs[maskPassX])/length(maskPassX)
            return(c(tprX,ppvX))
        }))

        shuffTPrPPVMAt <- matrix(unlist(listShuffTprPpvI),ncol=2,byrow=TRUE)

        return(shuffTPrPPVMAt)
    })

    shuffTprMat <- matrix(unlist(lapply(listShuffTprPPV,function(matX){
        return(matX[,1])
    })),ncol=nShuff)

    shuffPpvMat <- matrix(unlist(lapply(listShuffTprPPV,function(matX){
        return(matX[,2])
    })),ncol=nShuff)

    shuffTPrQs <- t(apply(shuffTprMat,1,function(rowX){
        return(quantile(rowX,c(0.01,0.5,0.99)))
    }))

    shuffPPvQs <- t(apply(shuffPpvMat,1,function(rowX){
        return(quantile(rowX,c(0.01,0.5,0.99)))
    }))

    colnames(shuffTPrQs) <- paste("shuffTPR",c("q01","q50","q99"),sep="_")
    colnames(shuffPPvQs) <- paste("shuffPPV",c("q01","q50","q99"),sep="_")

    #calculate area under the curve and difference in area from shuffles
    totalAUC <- calcAUC(xVals=actualTPrPPVMAt[,1],
        yVals=actualTPrPPVMAt[,2])
    shuffMedianAUC <- calcAUC(xVals=shuffTPrQs[,2],
        yVals=shuffPPvQs[,2])
    deltaAUC <- totalAUC - shuffMedianAUC
    normDeltaAUC <- deltaAUC/(1-shuffMedianAUC)

    allXs <- c(actualTPrPPVMAt[,1],unlist(shuffTPrQs))
    allYs <- c(actualTPrPPVMAt[,2],unlist(shuffPPvQs))

    if(makePlot){
        plot(x=allXs,y=allYs,type='n',xlab="TPR",ylab="apparent PPV")
        lines(x=actualTPrPPVMAt[,1],y=actualTPrPPVMAt[,2],
            type='l',col='black')

        seethruRed <- changeColorAlpha("red",newAlpha=100)
        polygon(x=c(shuffTPrQs[,1],rev(shuffTPrQs[,1])),
            y=c(shuffPPvQs[,1],rev(shuffPPvQs[,3])),
            col=seethruRed,border=NA)
        lines(x=shuffTPrQs[,2],y=shuffPPvQs[,2],col="red")
    }
    outList <- list()
    outList$actualTprPpvMat <- actualTPrPPVMAt
    outList$shuffTprQMat <- shuffTPrQs
    outList$shuffPpvQMat <- shuffPPvQs
    outList$AUC <- totalAUC
    outList$deltaAUC <- deltaAUC
    outList$normDeltaAUC <- normDeltaAUC

    invisible(outList)
}

#' Determine Optimal Distance Bin Size
#'
#' Analyzes the associations between DEG and CRE \link[S4Vectors]{Hits} to
#' determine the optimal distance bin size for further analysis.
#'
#' @param degCreHits A \link[S4Vectors]{Hits} object containing the
#' associations between DEG (Differentially Expressed Genes) and CRE
#' (Cis-Regulatory Element) hits with distances.
#' @param CreP A numeric vector of CRE p-values corresponding to the
#' associations in \code{degCreHits}.
#' @param fracMinKsMedianThresh Numeric value from 0 to 1 of the threshold
#' for minimum Kolmogorov-Smirnov (KS) median range. Determines the range of
#' KS statistics considered for optimal bin size. (Default: \code{0.2})
#' @param smallestTestBinSize Integer minimum number of associations in
#' each test bin. (Default: \code{100})
#' @param verbose Logical indicating whether to display progress messages.
#' (Default: \code{TRUE})
#'
#' @return A list containing:
#' \describe{
#'   \item{pickedBinSize}{The optimal distance bin size selected based on KS
#'   statistics.}
#'   \item{crePKsMat}{A matrix of distance bin sizes and their corresponding
#'   median KS statistics.}
#' }
#'
#' @details
#' Not exported. This function analyzes the associations between DEG and CRE
#' hits to determine the optimal distance bin size. It uses Kolmogorov-Smirnov
#' (KS) statistics to assess the difference in distribution between CRE
#' p-values for different distance bin sizes versus the global. The function
#' selects the largest bin size that falls within the specified fraction of
#' the KS median range.
#' This function operates within \link{runDegCre} on controlled inputs. It
#' will not run well on unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage:
#'
#' optimalBinSize <- distBinHeuristic(degCreHits = myDegCreHits,
#'                                    CreP = myCreP)
#'
#' # Access the selected bin size:
#' selectedBinSize <- optimalBinSize$pickedBinSize
#'
#' # Access the matrix of bin sizes and median KS statistics:
#' binStatsMatrix <- optimalBinSize$crePKsMat
#'
#' # Plot the results:
#' plot(binStatsMatrix[, 1], binStatsMatrix[, 2], type = "l", xlab =
#' "Distance Bin Size", ylab = "Median KS Statistic")
#' }
#'
#' @author Brian S. Roberts
#'
distBinHeuristic <- function(degCreHits,
                              CreP,
                              fracMinKsMedianThresh = 0.2,
                              smallestTestBinSize = 100,
                              verbose = TRUE) {

    sortPromCreQueryHits <-
        S4Vectors::queryHits(degCreHits)[order(S4Vectors::mcols(degCreHits)$assocDist)]

    sortPromCreSubjHits <-
        S4Vectors::subjectHits(degCreHits)[order(S4Vectors::mcols(degCreHits)$assocDist)]

    numAssocs <- length(sortPromCreSubjHits)

    if(verbose){
        message("Analyzing ",numAssocs,
            " DegGR to CreGR hits for optimal distance bin size.")
    }

    allSortHitIndices <- seq_len(numAssocs)

    allHitCrePs <- CreP[sortPromCreSubjHits]

    #the max number of test bins is set by input parameter "maxDistBins"
    #Bin size has units of number of hits

    minLowNumDistBins <- 2
    maxLowNumDistBins <- 100

    logLowTestNumDistBins <- seq(from=log10(minLowNumDistBins),
        to=log10(maxLowNumDistBins),
        length.out=20)

    highestNumDistBins <- floor(numAssocs/smallestTestBinSize)

    logHiTestNumDistBins <- seq(from=log10(2*maxLowNumDistBins),
        to=log10(0.5*highestNumDistBins),
        length.out=4)

    logHiTestNumDistBins <- c(logHiTestNumDistBins,log10(highestNumDistBins))

    logTestNumDistBins <- c(logLowTestNumDistBins,logHiTestNumDistBins)

    testNumDistBins <- unique(floor(10^logTestNumDistBins))

    testBinSizes <-
        unique(floor(numAssocs/testNumDistBins))

    crePKsMedianMat <- calcKStestStatMedian(testBinSizes=testBinSizes,
         allSortHitIndices=allSortHitIndices,
         allHitPs=allHitCrePs)

    outMedianKS <- crePKsMedianMat[,2]

    #pick the largest bin size that is within fracMinSSEThresh range of min
    KSMedianRange <- max(outMedianKS) - min(outMedianKS)

    pickedThresh <- min(outMedianKS) + fracMinKsMedianThresh*KSMedianRange

    maskPass <- which(outMedianKS<=pickedThresh)
    #cannot pick only 1 bin
    if(max(maskPass)==1){
        maskPass <- 2
    }

    pickedBinSize <-
        min(crePKsMedianMat[maskPass,1])

    colnames(crePKsMedianMat) <- c("distBinSize","medianBinKSstat")

    #return pickedBinSize and crePKsMedianMat for plotting
    outList <- list()
    outList$pickedBinSize <- pickedBinSize
    outList$crePKsMat <- crePKsMedianMat
    return(outList)
}

#' Calculate the Median KS Statistic by Distance Bin Size
#'
#' Calculates the median Kolmogorov-Smirnov (KS) statistic for different
#' distance bin sizes.
#'
#' @param testBinSizes A vector of integer values specifying the distance bin
#' sizes to be tested.
#' @param allSortHitIndices A vector of sorted hit indices.
#' @param allHitPs A numeric vector of p-values corresponding to the
#' associations.
#'
#' @return A matrix containing two columns:
#' \describe{
#'   \item{distBinSize}{The tested distance bin sizes.}
#'   \item{KSRMSE}{The median KS statistic for each distance bin size.}
#' }
#'
#' @details
#' Not exported. This function calculates the median Kolmogorov-Smirnov
#' (KS) statistic for different distance bin sizes. It splits the associations
#' into bins of specified sizes, calculates the KS statistic for each bin, and
#' returns the median KS statistic for each bin size. The KS statistic
#' measures the maximum difference between two cumulative distribution
#' functions, providing a measure of the difference between CRE p-value
#' distributions within each bin and the global distribution.
#' It is meant to run within \link{distBinHeuristic}. It will not run well
#' on unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' binSizes <- c(10000, 5000, 200)
#' hitIndices <- seq_len(1e5)
#' pValues <- runif(1e5)
#'
#' resultMatrix <- calcKStestStatMedian(testBinSizes = binSizes,
#'                                       allSortHitIndices = hitIndices,
#'                                       allHitPs = pValues)
#' }
#'
#' @author Brian S. Roberts
#'
calcKStestStatMedian <- function(testBinSizes,
                                 allSortHitIndices,
                                 allHitPs){

    # Pre-compute an ecdf for the reference (allHitPs)
    allHitPCumProb <- rank(allHitPs)/length(allHitPs)

    listKsMeanByBinSize <- lapply(testBinSizes, function(binSizeX){
        sortHitChunksListX <-
          split(allSortHitIndices,
                ceiling(seq_along(allSortHitIndices)/(binSizeX)))

        # Don't use the last bin for analysis because it is always smaller
        sortHitChunksListX <-
          sortHitChunksListX[seq_len((length(sortHitChunksListX)-1))]

        chunkKSstats <- unlist(lapply(sortHitChunksListX, function(hitIndicesY){
            hitYCrePs <- allHitPs[hitIndicesY]

            ksStatY <- fastKS(testSet=hitYCrePs,
                              testIndices=hitIndicesY,
                              refCumProbs=allHitPCumProb)
            return(ksStatY)
        }))

        kSMean <- median(chunkKSstats)

        return(kSMean)
    })

    KsMeanByBinSizeMat <- cbind(testBinSizes, unlist(listKsMeanByBinSize))
    colnames(KsMeanByBinSizeMat) <- c("distBinSize", "KSRMSE")
    return(KsMeanByBinSizeMat)
}

#' Perform Kolmogorov-Smirnov (KS) Test Without Calculating P-Value
#'
#' This function performs a Kolmogorov-Smirnov (KS) test quickly without
#' calculating the p-value. It measures the maximum difference between
#' cumulative probability distributions of a test set and a reference set.
#'
#' @param testSet Numeric vector representing the values of the test set.
#' @param testIndices Indices of elements in the test set to consider.
#' @param refCumProbs Numeric vector representing the pre-computed cumulative
#' probabilities of the reference set.
#'
#' @return Numeric value of the KS test statistic.
#'
#' @details
#' Not exported. The function compares the cumulative probability
#' distributions of the specified test set to the pre-computed cumulative
#' probabilities of the reference set. It returns the the KS statistic without
#' calculating the p-value (for computational speed).
#' It is meant to run within \link{calcKStestStatMedian}. It will not run well
#' on unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage of the function.
#' ks_stat <- fastKS(testSet = testData,
#'                             testIndices = indices,
#'                             refCumProbs = refCumProbs)
#' }
#'
#' @author Brian S. Roberts
#'
fastKS <- function(testSet, testIndices, refCumProbs) {
  testCumPs <- rank(testSet) / length(testSet)

  refCumPs <- refCumProbs[testIndices]

  testToRefDiffs <- abs(testCumPs - refCumPs)

  return(max(testToRefDiffs))
}

#' Get Associations and Distances Between Genomic Regions
#'
#' This function finds associations and distances between two sets of genomic
#' regions.
#'
#' @param DegGR A \link[GenomicRanges]{GRanges} object representing DEG TSSs.
#' @param CreGR A \link[GenomicRanges]{GRanges} object representing the CREs.
#' @param maxDist Integer value representing the maximum distance allowed for
#' associations. Regions further apart than this threshold will not be
#' associated. (Default: \code{1e6})
#'
#' @return A Hits object containing associations and their distances between
#' the genomic regions represented by DegGR and CreGR.
#'
#' @details
#' This function identifies associations between genomic regions from two
#' GenomicRanges objects (DegGR and CreGR) based on their spatial overlap
#' within a specified maximum distance threshold using
#' \link[GenomicRanges]{findOverlaps}. It calculates the distances with
#' \link[GenomicRanges]{distance} between associated regions and stores them
#' in the metadata of the \link[S4Vectors]{Hits} object.
#' Large values \code{maxDist} will require more computational resources.
#'
#' @examples
#' #Load sample data.
#' data(DexNR3C1)
#'
#' # Get hits with association distances.
#' hits <- getAssocDistHits(DegGR = DexNR3C1$DegGR,
#'                          CreGR = DexNR3C1$CreGR,
#'                          maxDist = 1e6)
#'
#' @author Brian S. Roberts
#'
#' @export
getAssocDistHits <- function(DegGR, CreGR, maxDist = 1e6) {
  hitsPromsToCRE <- IRanges::findOverlaps(DegGR, CreGR, maxgap = maxDist,
                                                ignore.strand = TRUE)

  # Get the association distances
  assocDist <-
    IRanges::distance(DegGR[S4Vectors::queryHits(hitsPromsToCRE)],
    CreGR[S4Vectors::subjectHits(hitsPromsToCRE)])

  distDf <- data.frame(assocDist = assocDist)

  S4Vectors::mcols(hitsPromsToCRE) <- distDf
  return(hitsPromsToCRE)
}

#' Calculate Enrichment Statistics for Dependent and Independent Data
#'
#' This function calculates enrichment statistics for dependent and
#' independent data based on the provided data frame and parameters.
#'
#' @param hitsWithDistDf A \link[S4Vectors]{DataFrame} derived from a Hits
#' object.
#' @param subHitsIndex Indices of the row subsets of \code{hitsWithDistDf} to
#' analyze.
#' @param dependPadj Numeric vector of adjusted p-values for the dependent
#' data.
#' @param independP Numeric vector of p-values for the independent data.
#' @param alpha Numeric significance level threshold for DEGs.
#'
#' @return A matrix containing calculated statistics for enrichment analysis,
#' including independent p-values, associated probabilities, and the total
#' number of observations.
#'
#' @details
#' Not exported. This function calculates enrichment statistics for dependent
#' and independent data based on provided adjusted p-values for dependent data
#' and p-values for independent data. It computes the associated probabilities,
#' which represent the probability of observing a significant association for
#' each set of data under the specified significance level threshold.
#' the independent variable is DegCre calculations is the CreP and the
#' dependent is the DEG adjusted p-values.
#' It is meant to run within \link{runDegCre}. It will not run well on
#' unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get stat results.
#' statsMatrix <- calcDependIndependEnrichStats(hitsWithDistDf = myHitsDf,
#'                                              subHitsIndex = mySubHits,
#'                                              dependPadj = myDependPadj,
#'                                              independP = myIndependP,
#'                                              alpha = 0.05)
#' }
#'
#' @author Brian S. Roberts
#'
calcDependIndependEnrichStats <- function(hitsWithDistDf, subHitsIndex,
                                          dependPadj,
                                          independP,
                                          alpha) {

  totalObs <- nrow(hitsWithDistDf[subHitsIndex,,drop=FALSE])

  subjHitsTemp <- hitsWithDistDf[subHitsIndex,2]
  queryHitsTemp <- hitsWithDistDf[subHitsIndex,1]

  # To run faster, run across approximately equal independent p-values

  log10P <- log10(independP[subjHitsTemp])
  log10RoundedP <- round(log10P, digits = 1)

  listMapLog10RoundedAllToUniq <- tapply(seq_along(log10RoundedP),
                                         INDEX = log10RoundedP,
                                         FUN = c)

  uniqlog10RoundedP <- as.numeric(names(listMapLog10RoundedAllToUniq))

  mapUniqToLog10RoundedAll <- unlist(listMapLog10RoundedAllToUniq)

  mapUniqRepTimes <- lapply(listMapLog10RoundedAllToUniq, length)

  # Filter to hit_df by p threshold
  AssocProbByPThreshUniq <- unlist(lapply(uniqlog10RoundedP,function(threshPX){
    maskHitsPass <- which(log10RoundedP <= threshPX)
    # The total number of "calls" is the length of maskHitsPass
    allPosCalls <- length(maskHitsPass)

    maskAssocQueries <- queryHitsTemp[maskHitsPass]
    numPassAlpha <- length(which(dependPadj[maskAssocQueries] <= alpha))

    truePosCalls <- numPassAlpha * (1 - alpha)

    expectProb <- truePosCalls / allPosCalls
    return(expectProb)
  }))

  # Convert AssocProbByPThreshUniq to the original (un-unique) set
  AssocProbByPThresh <- numeric(length = length(log10RoundedP))
  AssocProbByPThresh[mapUniqToLog10RoundedAll] <-
    rep(AssocProbByPThreshUniq, times = mapUniqRepTimes)

  allPMat <- cbind(independP[subjHitsTemp], AssocProbByPThresh, totalObs)

  colnames(allPMat) <- c("independP", "assocProb", "totalObs")

  return(allPMat)
}

#' Correct Association Probabilities
#'
#' This function corrects association probabilities based on distance bins and
#' reference association probabilities.
#'
#' @param sortHitsDf A \link[S4Vectors]{DataFrame} containing sorted hits data.
#' @param assocProbs A numeric vector of association probabilities.
#' @param refAssocProbs A numeric vector of reference association
#' probabilities. (Default: \code{NULL}, uses \code{assocProbs} if \code{NULL})
#'
#' @return A numeric vector of corrected association probabilities.
#'
#' @details
#' Not exported. This function corrects association probabilities within the
#' same distance bin based on reference association probabilities in lower
#' distance bins. It calculates adjusted association probabilities and
#' reference association probabilities for each distance bin and updates the
#' original association probabilities accordingly.
#' The principle is that for all associations involving a single CRE, those
#' associations to significant DEGs that span the shortest distances should be
#' weighted higher than those that span farther distances.
#' It is meant to run within \link{runDegCre}. It will not run well on
#' unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Distance bin correct association probabilities.
#' correctedProbs <- correctAssocProbs(sortHitsDf = mySortedHits,
#'                                     assocProbs = myAssocProbs,
#'                                     refAssocProbs = myRefAssocProbs)
#' }
#'
#' @author Brian S. Roberts
#'
correctAssocProbs <- function(sortHitsDf, assocProbs, refAssocProbs = NULL) {
  if (is.null(refAssocProbs)) {
    refAssocProbs <- assocProbs
  }

  # sortHitsDf[,4] is the dist bin Id. Assoc probs within the same distance
  # bin are only corrected by Assoc probs in lower dist bins.

  UniqSubjHitAdjProbs <- unlist(tapply(seq_along(assocProbs),
                                       INDEX = sortHitsDf[,2],
                                       FUN = function(indexX){

    assocProbsX <- assocProbs[indexX]
    refAssocProbsX <- refAssocProbs[indexX]

    assocDistBinsX <- sortHitsDf[indexX,4]
    uniqDistBinsX <- unique(assocDistBinsX)

    adjAssocProbsX <- assocProbsX
    adjRefAssocProbsX <- refAssocProbsX

    for (uniqDBinY in uniqDistBinsX) {
      subMaskLowerBinsY <- which(assocDistBinsX < uniqDBinY)
      subMaskBinY <- which(assocDistBinsX == uniqDBinY)

      adjAssocProbsY <- adjAssocProbsX[subMaskBinY]
      adjRefAssocProbsY <- adjRefAssocProbsX[subMaskBinY]

      if (length(subMaskLowerBinsY) > 0) {
        sumLowerProbs <- sum(adjRefAssocProbsX[subMaskLowerBinsY])
        adjAssocProbsY <- (adjAssocProbsY^2) / (sumLowerProbs + adjAssocProbsY)
        adjRefAssocProbsY <-
          (adjRefAssocProbsY^2) / (sumLowerProbs + adjRefAssocProbsY)

        maskNotFinite <- which(!is.finite(adjAssocProbsY))
        if (length(maskNotFinite) > 0) {
          adjAssocProbsY[maskNotFinite] <- 0
        }

        maskNotFiniteRef <- which(!is.finite(adjRefAssocProbsY))
        if (length(maskNotFiniteRef) > 0) {
          adjRefAssocProbsY[maskNotFiniteRef] <- 0
        }
      }

      adjAssocProbsX[subMaskBinY] <- adjAssocProbsY
      adjRefAssocProbsX[subMaskBinY] <- adjRefAssocProbsY
    }
    return(adjAssocProbsX)
  }))

  MapUniqSubjToOrig <- unlist(tapply(seq_along(assocProbs),
                                     sortHitsDf[,2],
                                     function(Ix){
    return(Ix)
  }))

  adjustProbs <- numeric(length = length(assocProbs))
  adjustProbs[MapUniqSubjToOrig] <- UniqSubjHitAdjProbs
  adjustProbs[which(assocProbs == 0)] <- 0

  return(adjustProbs)
}

#' Calculate Binomial FDR per Distance Bin
#'
#' Calculates the False Discovery Rate (FDR) for association probabilities
#' within a given distance bin using a binomial distribution approach.
#'
#' @param allDistBinsStatsMat A matrix containing statistics for all distance
#' bins.
#' @param chunkI An integer vector specifying the indices of the current
#' distance bin.
#' @param alphaVal Numeric from 0 to 1 of the DEG significance level.
#'
#' @return A numeric vector of FDR values for the specified distance bin.
#'
#' @details
#' Not exported. This function calculates the FDR for association probabilities
#' within a given distance bin. It uses a binomial distribution approach
#' (via \link[stats]{pbinom}) to estimate the FDR based on the number of
#' significant associations and the total number of associations in the bin.
#' Additionally, it adjusts FDR values for associations with low probabilities
#' or ties in the significance ranks.
#' It is meant to run within \link{runDegCre}. It will not run well on
#' unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get FDR from binomila distribtuion.
#' binomFDR <- calcBinomFDRperBin(allDistBinsStatsMat = myStatsMatrix,
#' chunkI = myChunkIndices, alphaVal = 0.05)
#' }
#'
#' @author Brian S. Roberts
#'

calcBinomFDRperBin <- function(allDistBinsStatsMat, chunkI, alphaVal) {

  chunkDegPadjX <- allDistBinsStatsMat[chunkI,5]
  chunkCrePX <- allDistBinsStatsMat[chunkI,3]

  chunkAdjAsscProbX <- allDistBinsStatsMat[chunkI,1]

  # Total expected DEGs is the number that pass alpha times (1-alpha)

  # Total DEGs
  totDegs <- length(which(chunkDegPadjX <= alphaVal)) * (1 - alphaVal)
  totAssocs <- length(chunkCrePX)

  successProb <- totDegs / totAssocs

  # CrepRanks is the number of balls drawn
  crepRanks <- rank(chunkCrePX, ties.method = "max")

  # DegHits is the number of white balls drawn
  degHits <- round(chunkAdjAsscProbX * crepRanks, 0)

  rankHitMat <- cbind(crepRanks, degHits)

  binomFDRX <- apply(rankHitMat, 1, function(rowY) {
    fdrY <- 1 - pbinom(q = rowY[2], size = rowY[1], prob = successProb)
    return(fdrY)
  })

  # The FDRs at the low end of the CreP ranks (most significant) must be
  # adjusted
  # because the number of trials in the binomial test is low for these.
  # The 10 most significant Cre Ps that have adjAsscProbs > successProb will
  # all get the FDR of the 10th (since it is powered)
  # However, if there are very few associations that pass FDR, then need
  # to pick a number lower than 10 to avoid problems

  maskAboveNull <- which(chunkAdjAsscProbX > successProb)

  if(length(maskAboveNull)==0){
    binomFDRX <- binomFDRX
  }
  else{
      if (length(maskAboveNull) > 50) {
        topFDRI <- 10
      } else {
        topFDRI <- ceiling(length(maskAboveNull) / 5)
        if (topFDRI == 0) {
          topFDRI <- 1
        }
      }
      nonNullCrepRanks <- crepRanks[maskAboveNull]
      rank10NonNull <- sort(nonNullCrepRanks)[topFDRI]

      maskInTop10 <- which(crepRanks <= rank10NonNull)

      maskInTop10NonNull <- intersect(maskAboveNull, maskInTop10)

      # Use the lowest FDR of any in the bottom 10
      lowerEndFDR <- min(binomFDRX[maskInTop10])

      binomFDRX[maskInTop10NonNull] <- lowerEndFDR
  }
  
  # Also, all adjAssocProbs of 0 should get FDR of 1
  #sometimes they don't because of weirdness with binom calc
  #this happens when many or all assocProb are 0
  maskZero <- which(chunkAdjAsscProbX == 0)
  
  binomFDRX[maskZero] <- 1
  
  return(binomFDRX)
}

#' Calculate Area Under the Curve (AUC)
#'
#' Calculates the Area Under the Curve (AUC) for a given set of x and y values
#' using the trapezoidal rule.
#'
#' @param xVals A numeric vector of x-values.
#' @param yVals A numeric vector of corresponding y-values.
#'
#' @return A numeric value representing the AUC.
#'
#' @details
#' Not exported. This function calculates the AUC for a given set of x and y
#' values using the trapezoidal rule. It provides a measure of the area under
#' the curve formed by the x and y values, which is often used to assess the
#' performance of models or the shape of a curve.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get AUC of quadratic curve.
#' x <- seq_len(10)
#' y <- x^2
#' auc <- calcAUC(x, y)
#' }
#' @author Brian S. Roberts
#'
calcAUC <- function(xVals, yVals) {
  auc <- abs(sum(diff(xVals) * (head(yVals, -1) + tail(yVals, -1)) / 2))
  return(auc)
}

