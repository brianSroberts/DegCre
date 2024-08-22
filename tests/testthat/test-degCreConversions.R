test_that("calcRawAssocProbOR", {
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

  calcORvec <- calcRawAssocProbOR(degCreResListDexNR3C1)

  # indices are based on top quantiles of assocProb
  testedIndices <- c(881,1003,15715,18745,19170,29471)
  
  calcORvec[testedIndices]
  
  testORs <- c(5.299773,6.270764,5.299773,2.725314,2.725314,4.245307)

  expect_equal(calcORvec[testedIndices],testORs,tolerance=1e-5)
})


test_that("convertdegCreResListToGInteraction", {
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

  calcGInter <-
    convertdegCreResListToGInteraction(degCreResList=degCreResListDexNR3C1)

  #pick subset of assocProb that span the range
  testGInterIndices <- c(40,59,10,22,373,555)
  mcols(calcGInter)$assocProb[testGInterIndices]
  
  testGInterAssocProbs <- c(0.1328839,0.5261069,0.1510457,0.3840299,
                            0.2136784,0.6242492)

  expect_equal(mcols(calcGInter)$assocProb[testGInterIndices],
               testGInterAssocProbs,
               tolerance=1e-5)

  #check mapping
  testHits <- degCreResListDexNR3C1$degCreHits
  maskGInterPassFDR <- which(mcols(testHits)$assocProbFDR<=0.05)

  passHits <- testHits[maskGInterPassFDR]

  passInDegGeneSymbs <- degCreResListDexNR3C1$DegGR$GeneSymb[queryHits(passHits)]

  passInCrePval <- degCreResListDexNR3C1$CreGR$pVal[subjectHits(passHits)]

  expect_identical(mcols(calcGInter)$Deg_GeneSymb,passInDegGeneSymbs)

  expect_equal(mcols(calcGInter)$Cre_pVal,passInCrePval)
})


test_that("convertDegCreDataFrame", {
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

  degCreDf <- convertDegCreDataFrame(degCreResList=degCreResListDexNR3C1,
                                     assocAlpha = 0.05)

  #pick subset of assocProb that span the range
  testIndices <- c(57,231,358,401,477)

  testAssocProbs <- c(0.6216533,0.1718222,0.5160128,0.3712580,0.4490510)

  expect_equal(degCreDf$assocProb[testIndices],
               testAssocProbs,
               tolerance=1e-5)

  #check mapping
  testHits <- degCreResListDexNR3C1$degCreHits
  maskPassFDR <- which(mcols(testHits)$assocProbFDR<=0.05)

  passHits <- testHits[maskPassFDR]

  passInDegGeneSymbs <- degCreResListDexNR3C1$DegGR$GeneSymb[queryHits(passHits)]

  passInCrePval <- degCreResListDexNR3C1$CreGR$pVal[subjectHits(passHits)]

  expect_identical(degCreDf$Deg_GeneSymb,passInDegGeneSymbs)

  expect_equal(degCreDf$Cre_pVal,passInCrePval)
})
