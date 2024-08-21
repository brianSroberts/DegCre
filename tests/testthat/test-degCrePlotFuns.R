test_that("plotDegCreBinHeuristics", {
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

  testPickBin <- 2277
  names(testPickBin) <- "distBinSize"
  expect_equal(calcPickBin,testPickBin)
})

test_that("plotDegCreAssocProbVsDist", {
  
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

  assocDistMat <-
    plotDegCreAssocProbVsDist(degCreResList=degCreResListDexNR3C1)

  dev.off()

  # Check if the plot file was created
  expect_true(file.exists(pdf_file),
              "Assoc vs dist plot file should be created")

  #Clean up the temporary plot file
  unlink(pdf_file)

  # check assocDistMat
  testRows <-c(1,3,5,10,12,16)
  testMedians <- c(0.4627001,0.2445692,0.2360954,0.1623336,0.1862349,
                   0.1997442)

  expect_equal(assocDistMat[testRows,3],testMedians,tolerance=1e-5)
})

test_that("getDistBinNullAssocProb", {
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

  # run getDistBinNullAssocProb
  calcNullMat <- getDistBinNullAssocProb(degCreResList = degCreResListDexNR3C1)

  testRows <- c(1,2,4,6,12,14,16)
  testNullProbs <- c(0.22434783,0.14565217,0.14304348,0.11913043,0.07956522,
                     0.08086957,0.10280105)

  expect_equal(calcNullMat[testRows,2],testNullProbs,tolerance=1e-5)
})
