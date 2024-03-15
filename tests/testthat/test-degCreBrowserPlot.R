test_that("plotBrowserDegCre", {
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

  # run plotDegCreAssocProbVsDist and check plot
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)

  calcBrowserOuts <- plotBrowserDegCre(degCreResList=degCreResListDexNR3C1,
                                   geneName="ERRFI1",
                                   geneNameColName="GeneSymb",
                                   CreSignalName="NR3C1")

  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file),
              "Assoc vs dist plot file should be created")

  # Clean up the temporary plot file
  unlink(pdf_file)

  # check outputs
  # Check CRE signal score converted to a dense GRanges for plotting
  testPlotRegionDf <- data.frame(chr="chr1",
                                 start=7302483,
                                 end=9293894)
  testPlotRegionGR <- GenomicRanges::makeGRangesFromDataFrame(testPlotRegionDf)

  expect_true(calcBrowserOuts$plotRegionGR==testPlotRegionGR,
              "plotRegionGR should be as expected")

  testSignalIndices <- c(98,92,3,21,44,90,5,55,37,82)
  testSignalScores <- c(-5.331180e-05,1.434569e-05,9.994705e-05,1.475679e-04,
                        3.488836e+00,1.753670e+00,2.638824e+00,
                        2.115627e+00,3.053984e+00,5.928976e-05)

  expect_equal(calcBrowserOuts$creSignalPlotGR$score[testSignalIndices],
               testSignalScores,tolerance=1e-4)

  # Check values in calcBrowserOuts$assocGinter
  testBrowserGInIndices <- c(1,3,6,9,10,11,17,19)
  testBrowserGInAssocProbs <- c(0.07096774,0.44448980,0.27720000,0.21300203,
                                0.09838509,0.39600000,0.35255396,0.30697674)

  expect_equal(mcols(calcBrowserOuts$assocGinter)$assocProb[testBrowserGInIndices],
               testBrowserGInAssocProbs,tolerance=1e-5)
})
