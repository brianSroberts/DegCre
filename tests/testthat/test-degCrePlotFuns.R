test_that("plotDegCreBinHeuristics", {

  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  # run plotDegCreBinHeuristic and check plot
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)

  calcPickBin <- plotDegCreBinHeuristic(degCreResList=degCreResListDexNR3C1)

  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file),
              "Assoc per DEG plot file should be created")

  # Clean up the temporary plot file
  unlink(pdf_file)

  testPickBin <- 3896
  names(testPickBin) <- "distBinSize"
  expect_equal(calcPickBin,testPickBin)
})

test_that("plotDegCreAssocProbVsDist", {

  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  # run plotDegCreAssocProbVsDist and check plot
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)

  assocDistMat <-
    plotDegCreAssocProbVsDist(degCreResList=degCreResListDexNR3C1)

  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file),
              "Assoc vs dist plot file should be created")

  #Clean up the temporary plot file
  unlink(pdf_file)

  # check assocDistMat
  testRows <-c(37,18,5,180,64,160,55,17,57,15)
  testMedians <- c(0.07221322,0.11525735,0.19435583,0.08991211,0.09245690,
                   0.06187500,0.08176129,0.09222222,0.06708419,0.10704815)

  expect_equal(assocDistMat[testRows,3],testMedians,tolerance=1e-5)
})

test_that("getDistBinNullAssocProb", {

  # bring in test degCre inputs
  data("DexNR3C1")

  degCreResListDexNR3C1 <- DegCre::runDegCre(DegGR=DexNR3C1$DegGR,
                                             DegP=DexNR3C1$DegGR$pVal,
                                             DegLfc=DexNR3C1$DegGR$logFC,
                                             CreGR=DexNR3C1$CreGR,
                                             CreP=DexNR3C1$CreGR$pVal,
                                             CreLfc=DexNR3C1$CreGR$logFC,
                                             verbose=FALSE)

  # run getDistBinNullAssocProb
  calcNullMat <- getDistBinNullAssocProb(degCreResList = degCreResListDexNR3C1)

  testRows <- c(74,158,91,183,23,143,149,3,90,130)
  testNullProbs <- c(0.05437885,0.04929671,0.04980493,0.04497690,0.05996920,
                     0.05310832,0.06301848,0.10392967,0.07521561,0.05691992)

  expect_equal(calcNullMat[testRows,2],testNullProbs,tolerance=1e-5)
})
