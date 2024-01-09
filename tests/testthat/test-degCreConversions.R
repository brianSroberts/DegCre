test_that("calcRawAssocProbOR", {
  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  calcORvec <- calcRawAssocProbOR(degCreResListDexNR3C1)

  # indices are based on top quantiles of assocProb
  testedIndices <- c(540187,12655,13784,34542,797,11790,288558,18858,160113,
                     107838)

  testORs <- c(1.0000000,0.9708563,1.0000000,1.3783743,1.2893615,1.8154621,
               2.6749891,3.2723824,4.7409410,17.0131004)

  expect_equal(calcORvec[testedIndices],testORs,tolerance=1e-5)
})


test_that("convertdegCreResListToGInteraction", {
  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  calcGInter <-
    convertdegCreResListToGInteraction(degCreResList=degCreResListDexNR3C1)

  #pick subset of assocProb that span the range
  testGInterIndices <- c(1592,1486,1600,1243,2831,663,2604,2429,851,1001,1773)

  testGInterAssocProbs <- c(0.04196739,0.06374245,0.06992565,0.07557252,
                            0.08150943,0.08985149,0.09997649,0.11765525,
                            0.14480226,0.20435780,0.99000000)

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

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  degCreDf <- convertDegCreDataFrame(degCreResList=degCreResListDexNR3C1,
                                     assocAlpha = 0.05)

  #pick subset of assocProb that span the range
  testIndices <- c(1592,1486,1600,1243,2831,663,2604,2429,851,1001,1773)

  testAssocProbs <- c(0.04196739,0.06374245,0.06992565,0.07557252,
                            0.08150943,0.08985149,0.09997649,0.11765525,
                            0.14480226,0.20435780,0.99000000)

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
