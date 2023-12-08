test_that("plotBrowserDegCre", {

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
                                 start=7264791,
                                 end=9293894)
  testPlotRegionGR <- makeGRangesFromDataFrame(testPlotRegionDf)

  expect_true(calcBrowserOuts$plotRegionGR==testPlotRegionGR,
              "plotRegionGR should be as expected")

  testSignalIndices <- c(98,92,3,21,44,90,5,55,37,82)
  testSignalScores <- c(-2.895593e-04,-8.009457e-05,3.268266e-05,5.022616e-05,
                        3.727415e+00,1.753670e+00,4.324037e-05,-9.424158e-05,
                        5.196334e+00,-3.871754e-05)

  expect_equal(calcBrowserOuts$creSignalPlotGR$score[testSignalIndices],
               testSignalScores,tolerance=1e-4)

  # Check values in calcBrowserOuts$assocGinter
  testBrowserGInIndices <- c(4,10,12,9,5,3,7,13,21,19)
  testBrowserGInAssocProbs <- c(0.42923077,0.19349146,0.11880000,0.13860000,
                                0.20915493,0.24750000,0.66000000,0.37125000,
                                0.06223495,0.16500000)

  expect_equal(mcols(calcBrowserOuts$assocGinter)$assocProb[testBrowserGInIndices],
               testBrowserGInAssocProbs,tolerance=1e-5)
})
