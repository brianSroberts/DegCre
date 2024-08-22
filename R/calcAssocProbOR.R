#' Calculate Odds-ratio
#'
#' Given a DegCre results list, this function calculates the odds-ratio of the 
#' association probability.
#'
#' @param degCreResList List of DegCre results.
#' @param type A character of either \code{"raw"} or \code{"adj"}.
#' (Default:\code{"adj"}).
#'
#' @return A numeric of the association probability odds-ratios
#'
#' @details
#' This function is similar to \link{calcRawAssocProbOR} and will mimic its 
#' function when \code{type = "raw"}. When \code{type = "raw"} the 
#' calculation operates on the \code{"rawAssocProb"} metadata.
#' When \code{type = "adj"} the calculation operates on the \code{"assocProb"} metadata. 
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
#' #Calculate odds-ratio.
#' 
#' calcOR <- calcAssocProbOR(degCreResListDexNR3C1)
#'
#' @author Brian S. Roberts
#'
#' @export

calcAssocProbOR <- function(degCreResList,
                            type="adj"){
  
  hitsX <- degCreResList$degCreHits
  alphaValX <- degCreResList$alphaVal
  degPadjX <- S4Vectors::mcols(hitsX)$DegPadj
  binIdsX <- S4Vectors::mcols(hitsX)$distBinId
  
  typeQuery <- "assocProb"
  if(type=="raw"){
    typeQuery <- "rawAssocProb"
  }
  
  maskType <- which(colnames(S4Vectors::mcols(hitsX)) == typeQuery)
  
  adjAssocProbsX <- S4Vectors::mcols(hitsX)[[maskType]]
  
  nullExpectedProbByBin <- tapply(degPadjX,binIdsX,function(degPadjY){
    expectTot <- (1-alphaValX)*length(which(degPadjY<=alphaValX))
    expectProbY <- expectTot/length(degPadjY)
    return(expectProbY)
  })
  uniqBinIds <- as.integer(names(nullExpectedProbByBin))
  nullExpectedProbByBin <- as.numeric(nullExpectedProbByBin)
  
  outAssocProbOR <- rep(NA,length(adjAssocProbsX))
  
  for(binIdY in uniqBinIds){
    maskY <- which(binIdsX == binIdY)
    maskBinY <- which(uniqBinIds == binIdY)
    
    nullExpectedProbBinY <- nullExpectedProbByBin[maskBinY]
    
    if(nullExpectedProbBinY == 0){
      outORY <- rep(1,length(maskY))
    }
    else{
      outORY <- adjAssocProbsX[maskY]/nullExpectedProbBinY
    }
    
    outAssocProbOR[maskY] <- outORY
  }
  return(outAssocProbOR)
}
