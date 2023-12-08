test_that("getExpectAssocPerDEG", {
  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  #run getExpectAssocPerDEG
  calcExpectAssocPerDegDf <-
    getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
                         geneNameColName = "GeneSymb")

  expect_equal(dim(calcExpectAssocPerDegDf),c(15644,11))

  # test at randomly sampled indices across non-trivial expectAssocs
  testedIndices <- c(313,91,198,5,422,52,195,176,354,269)
  testExpectAssocs <- c(1.0101684,3.6963609,2.0376538,11.4674811,0.1317811,
                        5.3510017,2.0764307,2.3116513,0.6852296,1.4027410)

  calcExpectAssocs <- calcExpectAssocPerDegDf$expectAssocs[testedIndices]

  expect_equal(calcExpectAssocs,testExpectAssocs,tolerance=1e-6)
})

test_that("plotExpectedAssocsPerDeg", {
  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  # run getExpectAssocPerDEG and check plot
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)

  calcExpectAssocPerDegDf <-
    getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
                         geneNameColName = "GeneSymb")

  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file),
              "Assoc per DEG plot file should be created")

  # Clean up the temporary plot file
  unlink(pdf_file)

  calcMedianExpectAssocs <- plotExpectedAssocsPerDeg(calcExpectAssocPerDegDf)

  expect_equal(calcMedianExpectAssocs,1.421911,tolerance=1e-6)
})
