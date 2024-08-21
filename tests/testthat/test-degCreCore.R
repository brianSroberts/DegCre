test_that("runDegCre", {
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- 
    DexNR3C1$DegGR[GenomicRanges::seqnames(DexNR3C1$DegGR)=="chr1"]
  subCreGR <- 
    DexNR3C1$CreGR[GenomicRanges::seqnames(DexNR3C1$CreGR)=="chr1"]

  # get runDegCre results with defaults.
  degCreResList <- DegCre::runDegCre(DegGR=subDegGR,
                             DegP=subDegGR$pVal,
                             DegLfc=subDegGR$logFC,
                             CreGR=subCreGR,
                             CreP=subCreGR$pVal,
                             CreLfc=subCreGR$logFC,
                             padjMethod = "qvalue")

  ## test results

  # DegCre list should have 5 slots
  expect_length(degCreResList,5)

  # DegCre Hits object should be 34169 hits long for this input
  expect_length(degCreResList$degCreHits,34169)

  # Check the picked distance bin size
  expect_equal(degCreResList$binHeurOutputs$pickedBinSize,2277)

  # Spot check selected association probs
  # These values were selected to test the dynamic range of non-zero
  # association probs and are thus the high end of the distribution
  # (90% and higher)

  testedIndices <- c(10613,18745,20768,25078,27163,30737)

  testedExpectedAssocProbs <- c(0.2109517,0.6114172,0.1306104,0.2327426,
                                0.4821006,0.1337257)

  calcAssocProbs <- mcols(degCreResList$degCreHits)$assocProb[testedIndices]

  expect_equal(calcAssocProbs,testedExpectedAssocProbs,tolerance=1e-4)

  # Spot check selected assocProbFDRs. These are selected from the same indices
  # as testedExpectedAssocProbs

  testedExpectedFDRs <- c(1.501954e-08,2.449965e-05,1.393260e-01,
                          1.080136e-11,0.000000e+00,1.000000e+00)

  calcFDRs <- mcols(degCreResList$degCreHits)$assocProbFDR[testedIndices]

  expect_equal(calcFDRs,testedExpectedFDRs,tolerance=1e-5)
})

test_that("optimizeAlphaDegCre", {
  #Checking the PR AUC vals used to pick optimal DEG alpha
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- 
    DexNR3C1$DegGR[GenomicRanges::seqnames(DexNR3C1$DegGR)=="chr1"]
  subCreGR <- 
    DexNR3C1$CreGR[GenomicRanges::seqnames(DexNR3C1$CreGR)=="chr1"]

  alphaTestSet <- c(0.001,0.003,0.005,0.01,0.05,0.1)

  calcAlphaOptList <- optimizeAlphaDegCre(DegGR=subDegGR,
                                          DegP=subDegGR$pVal,
                                          DegLfc=subDegGR$logFC,
                                          CreGR=subCreGR,
                                          CreP=subCreGR$pVal,
                                          CreLfc=subCreGR$logFC,
                                          testedAlphaVals=alphaTestSet,
                                          padjMethod="bonferroni",
                                          verbose=FALSE)
  
  

  testAlphaAUCs <- c(0.01915195,0.02240950,0.02226499,0.02358739,
                     0.02764595,0.02628054)

  calcAUC <- calcAlphaOptList$alphaPRMat[,2]
  calcAUC <- unname(calcAUC)

  expect_equal(calcAUC,testAlphaAUCs,tolerance=1e-5)
})


test_that("degCrePRAUC", {
  # Checking invisible return of degCrePRAUC()
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

  expect_equal(testPrAUCList$AUC,0.07262013,tolerance=1e-6)

  # must have low tolreance for normDeltaAUC because it is compared to
  # a null distribution created by random sampling
  expect_equal(testPrAUCList$normDeltaAUC,0.06503126,tolerance=1e-2)
})
