test_that("plotDegCreBinHeuristics", {

  library(GenomicRanges)
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- DexNR3C1$DegGR[which(seqnames(DexNR3C1$DegGR)=="chr1")]
  subCreGR <- DexNR3C1$CreGR[which(seqnames(DexNR3C1$CreGR)=="chr1")]
  
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

  library(GenomicRanges)
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- DexNR3C1$DegGR[which(seqnames(DexNR3C1$DegGR)=="chr1")]
  subCreGR <- DexNR3C1$CreGR[which(seqnames(DexNR3C1$CreGR)=="chr1")]
  
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
  testMedians <- c(0.27720000,0.14430000,0.10913043,0.10129344,0.07563284,
                   0.08984293)

  expect_equal(assocDistMat[testRows,3],testMedians,tolerance=1e-5)
})

test_that("getDistBinNullAssocProb", {

  library(GenomicRanges)
  # bring in test degCre inputs
  data("DexNR3C1")
  
  subDegGR <- DexNR3C1$DegGR[which(seqnames(DexNR3C1$DegGR)=="chr1")]
  subCreGR <- DexNR3C1$CreGR[which(seqnames(DexNR3C1$CreGR)=="chr1")]
  
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
  testNullProbs <- c(0.18956522,0.12826087,0.12695652,0.10565217,0.06739130,
                     0.07521739,0.08984293)

  expect_equal(calcNullMat[testRows,2],testNullProbs,tolerance=1e-5)
})
