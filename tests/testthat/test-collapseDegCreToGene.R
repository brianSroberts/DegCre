test_that("collapseDegCreToGene", {
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- 
    DexNR3C1$DegGR[GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1"]
  subCreGR <- 
    DexNR3C1$CreGR[GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1"]
  
  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=subDegGR,
                                             DegP=subDegGR$pVal,
                                             DegLfc=subDegGR$logFC,
                                             CreGR=subCreGR,
                                             CreP=subCreGR$pVal,
                                             CreLfc=subCreGR$logFC,
                                             verbose=FALSE)
  
  collapseDegCreResList <- collapseDegCreToGene(degCreResListDexNR3C1)
  
  collapseHits  <- collapseDegCreResList$degCreHits
  
  testIndices <- c(1,100,1000,5000,10000,15000,17500,19989)
  
  expectVals <- c(1,51,222,711,1651,2321,2508,2784)
  
  expect_equal(queryHits(collapseHits)[testIndices],expectVals)
  
})