test_that("getExpectAssocPerDEG", {
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

  #run getExpectAssocPerDEG
  calcExpectAssocPerDegDf <-
    getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
                         geneNameColName = "GeneSymb")

  expect_equal(dim(calcExpectAssocPerDegDf),c(1616,11))

  # test at randomly sampled indices across non-trivial expectAssocs
  testedIndices <- c(1,3,9,12,28,35,59)
  testExpectAssocs <- c(41.240723,29.320480,15.459290,14.803924,8.900007,
                        7.281138,5.762932)

  calcExpectAssocs <- calcExpectAssocPerDegDf$expectAssocs[testedIndices]

  expect_equal(calcExpectAssocs,testExpectAssocs,tolerance=1e-4)
})

test_that("plotExpectedAssocsPerDeg", {
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

  # run getExpectAssocPerDEG and check plot
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)

  calcExpectAssocPerDegDf <-
    getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
                         geneNameColName = "GeneSymb")
  
  medianExpAssocs <- plotExpectedAssocsPerDeg(calcExpectAssocPerDegDf)
  
  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file),
              "Assoc per DEG plot file should be created")

  # Clean up the temporary plot file
  unlink(pdf_file)

  calcMedianExpectAssocs <- medianExpAssocs

  expect_equal(calcMedianExpectAssocs,2.414807,tolerance=1e-4)
})
