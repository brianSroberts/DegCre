test_that("convDegCreResListToCreGeneScoreGR", {
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
  
  creToGeneGR <- 
    DegCre::convDegCreResListToCreGeneScoreGR(degCreResListDexNR3C1,
                                              scoreType="assocProb",
                                              geneColname="GeneSymb",
                                              onlyDEGs=TRUE,
                                              DEgAlpha=NULL,
                                              degPadjColname="pAdj")
  
  testIndices <- c(1,2,5,100,500,1000,2000,3000,3500,3757)
  
  expectVals <- c(0.6180612,0.0951286,0.0263663,0.2223335,0.5040768,0.3558702,
                  0.0264783,0.3095202,0.3320311,0.2194128)
  
  expectGeneNames <- c("CHD5","CHD5","CHD5","ERRFI1","TRIM63","NFIA","MCL1",
                       "TMCC2","LEFTY1","GREM2")
  
  expect_equal(creToGeneGR$predictScore[testIndices],expectVals,tolerance=1e-4)
  
  expect_equal(creToGeneGR$predictGene[testIndices],expectGeneNames)
})
