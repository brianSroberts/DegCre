#' Collapse DegCre associations to a single TSS per gene
#'
#' Given a DegCre results list, this function finds associations between
#' the same CRE and multiple TSSs of the same gene and keeps the nearest
#' TSS only.
#'
#' @param degCreResList List of DegCre results.
#' @param method Method for choosing between multiple TSS. Currently only
#' supported for "nearest". (Default:\code{"nearest"})
#' @param geneColname The name of the metadata column in DegGR within
#' \code{degCreResList} that has the gene name (Default:\code{"GeneSymb"})
#' @param useParallel Logical whether to implement parallelization with
#' \link[BiocParallel]{BiocParallel} package. (Default:\code{FALSE})
#'
#' @return A \code{degCreResList} with only the shortest TSS to CRE
#' association per gene.
#'
#' @details
#' Often, the DegGR input to DegCre will contain multiple TSS's for a given
#' gene. DegCre will create associations to all of them. Often, downstream
#' analyses need to only have on CRE to TSS association per gene. This function
#' modifies the DegCre \link[S4Vectors]{Hits} to keep the shortest (minimum)
#' genomic distance association.
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
#' #Convert to single TSS per association.
#' 
#' degCreResListUniqTSS <- collapseDegCreToGene(degCreResListDexNR3C1,
#'                                              method = "nearest",
#'                                              geneColname = "GeneSymb",
#'                                              useParallel=FALSE)
#'
#' @author Brian S. Roberts
#'
#' @export
collapseDegCreToGene <- function(degCreResList,
                                 method = "nearest",
                                 geneColname = "GeneSymb",
                                 useParallel=TRUE) {
  
  BpParam <- BiocParallel::SerialParam()
  if(useParallel){
    #check for registered BiocParallel backend
    currentBpParam <- BiocParallel::bpparam()
    if(!inherits(currentBpParam,"SerialParam")){
      BpParam <- currentBpParam
    }else{
      BpParam <- BiocParallel::MulticoreParam()
    }
  }
  
  hitsX <- degCreResList$degCreHits
  
  maskGeneColName <- which(colnames(S4Vectors::mcols(degCreResList$DegGR)) ==
                             geneColname)
  
  hitsOrdGeneNames <- 
    S4Vectors::mcols(degCreResList$DegGR)[,maskGeneColName][queryHits(hitsX)]
  geneNamesSubjHitsHash <- paste0(hitsOrdGeneNames, "_", subjectHits(hitsX))
  
  # Compute the frequency of each hash
  rleGeneNameSubjHash <- rle(geneNamesSubjHitsHash)
  hashesNeedAttn <- rleGeneNameSubjHash$values[rleGeneNameSubjHash$lengths > 1]
  
  # Identify indices that need attention
  maskHashesNeedAttn <- which(geneNamesSubjHitsHash %in% hashesNeedAttn)
  
  if (length(maskHashesNeedAttn) < 1) {
    allKeepHitIndices <- seq_along(geneNamesSubjHitsHash)
    message("No DegCre associations found to collapse on gene.")
  } else {
    message("Processing associations for gene collapse")
    
    maskHashesDontNeedAttn <- 
      setdiff(seq_along(geneNamesSubjHitsHash), maskHashesNeedAttn)
    
    # Handle attention-needed hashes
    if (method == "nearest") {
      
      uniqGeneNamesSubjHitsHashNeedAttn <- 
        unique(geneNamesSubjHitsHash[maskHashesNeedAttn])
      
      attnHashesKeepIndices <- 
        BiocParallel::bplapply(uniqGeneNamesSubjHitsHashNeedAttn,
                               FUN = function(geneHashX) {
                                 hitsIndicesX <- 
                                   maskHashesNeedAttn[which(geneNamesSubjHitsHash[maskHashesNeedAttn] == geneHashX)]
                                 return(hitsIndicesX[which.min(S4Vectors::mcols(hitsX)$assocDist[hitsIndicesX])])
                               },
                               BPPARAM = BpParam
        )
      
      attnHashesKeepIndices <- unlist(attnHashesKeepIndices, use.names = FALSE)
    } else {
      stop("Unsupported method: ", method)
    }
    
    allKeepHitIndices <- sort(c(maskHashesDontNeedAttn, attnHashesKeepIndices))
  }
  
  degCreResList$degCreHits <- degCreResList$degCreHits[allKeepHitIndices]
  return(degCreResList)
}