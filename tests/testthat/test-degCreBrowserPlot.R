test_that("plotBrowserDegCre", {
  # bring in test degCre inputs
  set.seed(1976)
  
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
                                 start=6384796,
                                 end=9293894)
  testPlotRegionGR <- GenomicRanges::makeGRangesFromDataFrame(testPlotRegionDf)

  expect_true(calcBrowserOuts$plotRegionGR==testPlotRegionGR,
              "plotRegionGR should be as expected")

  testSignalIndices <- c(15,27,31,32,35,44,55,56,57,65)
  testSignalScores <- c(0.3323163,1.9551651,0.6089527,1.7602246,2.6388244,
                        1.0231081,0.0000000,5.1963341,3.0539838,6.2593598)

  expect_equal(calcBrowserOuts$creSignalPlotGR$score[testSignalIndices],
               testSignalScores,tolerance=1e-4)

  # Check values in calcBrowserOuts$assocGinter
  testBrowserGInIndices <- c(1,3,6,9,10,11,17,19,34,35)
  testBrowserGInAssocProbs <- c(0.5060245,0.1348534,0.6539693,0.1096556,
                                0.1489473,0.1444584,0.2487383,0.3334610,
                                0.5085380,0.4763748)
  
  actualVals <- 
    calcBrowserOuts$assocGinter$assocProb[testBrowserGInIndices]

  expect_equal(actualVals,testBrowserGInAssocProbs,tolerance=1e-4)
})
