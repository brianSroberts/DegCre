#' Convert a degCreResList to a creGeneScoreGR
#'
#' Converts a degCreResList to a hit indexed \link[GenomicRanges]{GRanges} with 
#' metadata columns of predicted gene target ("predictGene") and the association 
#' score ("score").
#'
#' @param degCreResList A \code{degCreResList}
#' @param scoreType Character of one of \code{assocProb}, \code{adjOR}, 
#' \code{rawAssocProb}. (Default:\code{assocProb})
#' See Details for a description of these options.
#' @param geneColname The name of the metadata column in \code{DegGR} in 
#' \code{degCreResList} that contains the gene names. 
#' (Default:\code{"GeneSymb"})
#' @param onlyDEGs Logical. If \code{TRUE}, only those associations involving
#' a gene with an adjusted p-value less than or equal to \code{DEgAlpha} are
#' returned. (Default:\code{TRUE})
#' @param DEgAlpha Significance threshold for DEGs in associations to report.
#' (Default:\code{NULL}.)
#' @param degPadjColname The metadata column name in \code{DegGR} that contains
#' the adjusted p-values. (Default:\code{pAdj})
#' 
#' @details This function simplifies a \code{degResList} to a 
#' \link[GenomicRanges]{GRanges} with metadata of the predicted associated gene
#' and a score of that prediction (\code{predictScore}). 
#' If \code{scoreType = "assocProb"}, the score is the DegCre association
#' probability. By setting \code{scoreType = "adjOR"} the odds-ratio of the 
#' association probability is calculated with \link{calcAssocProbOR} and
#' reported. This can be useful to emphasize associations
#' that are greatly above background and will tend to weight longer distance
#' associations more heavily. By setting \code{scoreType = "rawAssocProb"} the 
#' raw association probability is used. This is the probability of association
#' without considering the probability that the involved DEG is truly 
#' differential. This score may be useful when comparing to CRISPRi validation
#' in the absence of a perturbation.
#' 
#' @return A \link[GenomicRanges]{GRanges} of CRE regions with metadata columns
#' of \code{predictGene} and \code{predictScore}.
#' @author Brian S. Roberts
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
#' #Convert to GR with assocProb as score.
#' 
#' creToGeneGR <- convDegCreResListToCreGeneScoreGR(degCreResListDexNR3C1,
#'                                                  scoreType="assocProb",
#'                                                  geneColname="GeneSymb",
#'                                                  onlyDEGs=TRUE,
#'                                                  DEgAlpha=NULL,
#'                                                  degPadjColname="pAdj")
#'
#' @export
#' 
convDegCreResListToCreGeneScoreGR <- function(degCreResList,
                                              scoreType="assocProb",
                                              geneColname="GeneSymb",
                                              onlyDEGs=TRUE,
                                              DEgAlpha=NULL,
                                              degPadjColname="pAdj"){
  
  #get hits-indexed CreGR
  creHitsIndices <- S4Vectors::subjectHits(degCreResList$degCreHits)
  hitIndexedCreGR <- degCreResList$CreGR[creHitsIndices]
  
  hitIndexedMcolsDf <- 
    as.data.frame(S4Vectors::mcols(degCreResList$degCreHits))
  
  if(scoreType=="assocProb"){
    hitIndexedScore <- hitIndexedMcolsDf$assocProb
  }
  
  if(scoreType=="adjOR"){
    adjORs <- calcAssocProbOR(degCreResList)
    hitIndexedScore <- adjORs
  }
  
  if(scoreType=="rawAssocProb"){
    hitIndexedScore <- hitIndexedMcolsDf$rawAssocProb
  }
  
  degHitsIndices <- S4Vectors::queryHits(degCreResList$degCreHits)
  hitIndexedDegGR <- degCreResList$DegGR[degHitsIndices]
  
  mcolsDegHitIndexedDf <- as.data.frame(S4Vectors::mcols(hitIndexedDegGR))
  
  maskGeneCol <- which(colnames(mcolsDegHitIndexedDf)==geneColname)
  
  hitIndexedGeneName <- as.character(mcolsDegHitIndexedDf[,maskGeneCol])
  
  outCreGeneScoreGR <- granges(hitIndexedCreGR)
  
  outCreGeneScoreGR$predictGene <- hitIndexedGeneName
  outCreGeneScoreGR$predictScore <- hitIndexedScore
  
  if(onlyDEGs){
    if(is.null(DEgAlpha)){
      degAlpha <- degCreResList$alphaVal
    }else{
      degAlpha <- DEgAlpha
    }
    
    maskPadj <- which(colnames(mcolsDegHitIndexedDf) == degPadjColname)
    maskOnlyDegs <- which(mcolsDegHitIndexedDf[,maskPadj] <= degAlpha)
    outCreGeneScoreGR <- outCreGeneScoreGR[maskOnlyDegs]
  }
  
  return(outCreGeneScoreGR)
}