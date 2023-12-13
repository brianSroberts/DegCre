test_that("runDegCre", {
  # bring in test input data
  data("DexNR3C1")

  # get runDegCre results with defaults.
  degCreResList <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                             DegP=DexNR3C1$DegGR$pVal,
                             DegLfc=DexNR3C1$DegGR$logFC,
                             CreGR=DexNR3C1$CreGR,
                             CreP=DexNR3C1$CreGR$pVal,
                             CreLfc=DexNR3C1$CreGR$logFC)

  ## test results

  # DegCre list should have 5 slots
  expect_length(degCreResList,5)

  # DegCre Hits object should be 779231 hits long for this input
  expect_length(degCreResList$degCreHits,779231)

  # Check the picked distance bin size
  expect_equal(degCreResList$binHeurOutputs$pickedBinSize,3896)

  # Spot check selected association probs
  # These values were selected to test the dynamic range of non-zero
  # association probs and are thus the high end of the distribution
  # (95% and higher)

  testedIndices <- c(540187,12655,13784,34542,797,11790,288558,18858,160113,
                     107838)

  testedExpectedAssocProbs <- c(0.01654992,0.04539301,0.10392967,0.13274637,
                                0.16119667,0.22697005,0.29432432,0.35007605,
                                0.41441860,0.99000000)

  calcAssocProbs <- mcols(degCreResList$degCreHits)$assocProb[testedIndices]

  expect_equal(calcAssocProbs,testedExpectedAssocProbs,tolerance=1e-6)

  # Spot check selected assocProbFDRs. These are selected from the same indices
  # as testedExpectedAssocProbs

  testedExpectedFDRs <-
    mcols(degCreResList$degCreHits)$assocProbFDR[testedIndices]

  calcFDRs <- mcols(degCreResList$degCreHits)$assocProbFDR[testedIndices]

  expect_equal(calcFDRs,testedExpectedFDRs,tolerance=1e-6)
})

test_that("optimizeAlphaDegCre", {
  #Checking the PR AUC vals used to pick optimal DEG alpha

  # bring in test degCre inputs
  data("DexNR3C1")

  alphaTestSet <- c(0.0005,0.001,0.003,0.005,0.01,0.05,0.1,0.2)

  calcAlphaOptList <- optimizeAlphaDegCre(DegGR=DexNR3C1$DegGR,
                                          DegP=DexNR3C1$DegGR$pVal,
                                          DegLfc=DexNR3C1$DegGR$logFC,
                                          CreGR=DexNR3C1$CreGR,
                                          CreP=DexNR3C1$CreGR$pVal,
                                          CreLfc=DexNR3C1$CreGR$logFC,
                                          testedAlphaVals=alphaTestSet,
                                          verbose=FALSE)

  testAlphaAUCs <- c(0.008622081,0.008633184,0.009065857,0.009187438,
                     0.009699195,0.010220304,0.009113340,0.007995100)

  calcAUC <- calcAlphaOptList$alphaPRMat[,2]
  calcAUC <- unname(calcAUC)

  expect_equal(calcAUC,testAlphaAUCs,tolerance=1e-5)
})


test_that("degCrePRAUC", {
  # Checking invisible return of degCrePRAUC()

  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  # Capture the plot produced by the function
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)

  testPrAUCList <- degCrePRAUC(degCreResList=degCreResListDexNR3C1,
                               makePlot=TRUE)

  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file), "PR plot file should be created")

  # Clean up the temporary plot file
  unlink(pdf_file)

  ## testPrAUCList should have 6 slots
  expect_length(testPrAUCList,6)

  expect_equal(dim(testPrAUCList$actualTprPpvMat),c(200,2))
  expect_equal(dim(testPrAUCList$shuffTprQMat),c(200,3))
  expect_equal(dim(testPrAUCList$shuffPpvQMat),c(200,3))

  expect_equal(testPrAUCList$AUC,0.009699195,tolerance=1e-6)

  # must have low tolreance for normDeltaAUC because it is compared to
  # a null distribution created by random sampling
  expect_equal(testPrAUCList$normDeltaAUC,0.006629568,tolerance=1e-2)
})
