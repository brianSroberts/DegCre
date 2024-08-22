test_that("collapseDegCreToGene", {
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- 
    DexNR3C1$DegGR[GenomicRanges::seqnames(DexNR3C1$DegGR)=="chr1"]
  subCreGR <- 
    DexNR3C1$CreGR[GenomicRanges::seqnames(DexNR3C1$CreGR)=="chr1"]
  
  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=subDegGR,
                                             DegP=subDegGR$pVal,
                                             DegLfc=subDegGR$logFC,
                                             CreGR=subCreGR,
                                             CreP=subCreGR$pVal,
                                             CreLfc=subCreGR$logFC,
                                             verbose=FALSE)
  
  collapseDegCreResList <- collapseDegCreToGene(degCreResListDexNR3C1,
                                                useParallel=FALSE)
  
  collapseHits  <- collapseDegCreResList$degCreHits
  
  testIndices <- c(1,100,1000,5000,10000,15000,25000,34000)
  
  expectVals <- c(1,30,156,404,868,1573,2260,2775)
  
  expect_equal(queryHits(collapseHits)[testIndices],expectVals)
  
})