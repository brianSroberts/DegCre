test_that("calcAssocProbOR", {
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
  
  #with type="adj"
  assocProbOR <- calcAssocProbOR(degCreResListDexNR3C1,type="adj")

  testIndices <- c(1,2,5,111,112,116,201,249,280,3757)
  
  expectVals <- c(0.0000000,1.2664813,1.6039063,0.0000000,0.3280634,
                  0.9637354,0.2757159,2.7549238,0.2168754,2.7789968)
  
  expect_equal(assocProbOR[testIndices],expectVals,tolerance=1e-4)
  
  #with type="adj"
  rawAssocProbOR <- calcAssocProbOR(degCreResListDexNR3C1,type="raw")
  
  testIndices <- c(1,2,5,111,112,116,201,249,280,3757)
  
  rawExpectVals <- c(2.923578,2.923578,2.666548,2.991871,1.832532,2.587552,
                     1.942905,2.755082,2.084968,2.778999)
  
  expect_equal(rawAssocProbOR[testIndices],rawExpectVals,tolerance=1e-4)
})


