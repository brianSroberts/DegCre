test_that("calcRawAssocProbOR", {
  library(GenomicRanges)
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- DexNR3C1$DegGR[which(seqnames(DexNR3C1$DegGR)=="chr1")]
  subCreGR <- DexNR3C1$CreGR[which(seqnames(DexNR3C1$CreGR)=="chr1")]

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

  testORs <- c(3.044118,3.672581,3.044118,2.611239,2.611239,3.373333)

  expect_equal(calcORvec[testedIndices],testORs,tolerance=1e-5)
})


test_that("convertdegCreResListToGInteraction", {
  library(GenomicRanges)
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- DexNR3C1$DegGR[which(seqnames(DexNR3C1$DegGR)=="chr1")]
  subCreGR <- DexNR3C1$CreGR[which(seqnames(DexNR3C1$CreGR)=="chr1")]
  
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

  testGInterAssocProbs <- c(0.2368760,0.3641969,0.3611066,0.1606195,
                            0.4950000,0.2842414)

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
  library(GenomicRanges)
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- DexNR3C1$DegGR[which(seqnames(DexNR3C1$DegGR)=="chr1")]
  subCreGR <- DexNR3C1$CreGR[which(seqnames(DexNR3C1$CreGR)=="chr1")]

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

  testAssocProbs <- c(0.2368760,0.2842414,0.1607254,0.4950000,0.3641969)

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
