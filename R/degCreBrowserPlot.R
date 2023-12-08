#' Make Browser Plots from DegCre Results
#'
#' Creates browser plots of specified genomic regions or gene regions based on
#' the provided DegCre analysis results.
#'
#' @param degCreResList List of DegCre results.
#' @param assocAlpha Numeric value from 0 to 1 of significance threshold for
#' associations. (Default: \code{0.05})
#' @param browserWinPad Numeric value of the padding size (in base pairs) to
#' extend the plotting region. (Default: \code{1000})
#' @param geneName Character of name of the gene of interest. If specified,
#' the function will plot the region associated with this gene.
#' (Default: \code{NULL})
#' @param plotRegionGR \link[GenomicRanges]{GRanges} of length 1 specifying the
#' region to plot. If provided, geneName is ignored. (Default: \code{NULL})
#' @param CreSignalName Character of name of the differential CRE signal track.
#' For plot labeling purposes only (Default: \code{CRE})
#' @param assembly Character of genome assembly name, e.g., "hg38". Must be one
#' of accepted inputs to \code{assembly} argument to
#' \link[plotgardener]{plotGenomeLabel}. (Default: \code{hg38})
#' @param plotWidth Numeric value of width of the browser plot in inches.
#' (Default: \code{dev.size("in")[1]})
#' @param plotHeight  Numeric value of height of the browser plot in inches.
#' (Default: \code{dev.size("in")[2]})
#' @param plotXbegin Numeric value of the width of the left margin (where
#' track annotations will be) in inches. (Default: \code{0.9})
#' @param mergeGenePromotersDist Maximum distance (in base pairs) for merging
#' promoters of the same gene in plot. (Default: \code{1000})
#' @param sigPlotMaxY Numeric value of maximum value for the CRE differential
#' signal plot (Y-axis). (Default: \code{4})
#' @param assocColorRange Numeric vector of values from 0 to 1 of length 2.
#' These values specify the lower and upper values of DegCre association
#' probabilities for color saturation for arch color. (Default: \code{NULL}.
#' If \code{NULL} will be set to 0 and maximum association probability in
#' input data.)
#' @param lowAssocColor Character color for low saturation point of association
#' probabilities in arches plot. (Default: \code{#88CCEE})
#' @param hiAssocColor Character color for high saturation point of association
#' probabilities in arches plot. (Default: \code{#CC6677})
#' @param signalColor Character color for the CRE differential signal plot.
#' (Default: \code{#DDCC77})
#' @param geneLabelFontsize Numeric of font size (as implemented in
#' \href{\url{https://bioconductor.org/packages/release/bioc/html/plotgardener.html}}{plotgardener})
#' for gene labels. (Default: \code{8})
#' @param axisFontSize Numeric of font size for axis labels and tick marks.
#' (Default: \code{6})
#' @param panelTitleFontSize Numeric of font size for panel titles.
#' (Default: \code{7})
#' @param geneNameColName Character of name of the column in DegGR metadata
#' that was inputted to \link{runDegCre} that contains gene names to query by
#' \code{geneName}. (Default: \code{NULL}. If omitted, the column name will be
#' guessed, with warnings if not.)
#' @param geneHighlightDf \link[S4Vectors]{DataFrame} specifying genes to
#' highlight in the plot as accepted by \link[plotgardener]{plotGenes} argument
#' \code{geneHighlights}.
#' @param dePrioritizeSmallRNA Logical, indicating whether small RNA genes
#' should be deprioritized in plotting. (Default: \code{TRUE})
#' @param useLogFC Logical, indicating whether to use log-fold change values
#' for the CRE differential signal. (Default: \code{TRUE})
#' @param creSignalBinRes Bin resolution in base pairs for the CRE signal
#' track. Only used for initial calculation and will likely differ from display
#' resolution. (Default: \code{100})
#'
#' @importFrom plotgardener plotPairsArches annoHeatmapLegend colorby
#' plotSignal plotGenomeLabel plotGenes plotText pageCreate
#'
#' @details
#' This function uses
#' \href{\url{https://bioconductor.org/packages/release/bioc/html/plotgardener.html}}{plotgardener}
#' functionality to generate browser plots for visualizing DegCre analysis
#' results in specified regions. The user can input genomic regions or gene
#' names.
#' The output plot consists of an arches plot made with
#' \link[plotgardener]{plotPairsArches} of DegCre associations colored by
#' association probability.
#' Below is a signal track plot made via \link[plotgardener]{plotSignal} of the
#' data in \code{CreGR} that was inputted to \link{runDegCre}.
#' This plot displays the signed negative log p-value, meaning the
#' \eqn{-log_{10}(p_{CRE})} multiplied by the sign of the the CRE log
#' fold-change.
#' Beneath this panel genomic coordinates via
#' \link[plotgardener]{plotGenomeLabel} and gene models via
#' \link[plotgardener]{plotGenes} are displayed.
#'
#' @return Invisibly, a named list containing:
#' \describe{
#'   \item{plotRegionGR}{\link[GenomicRanges]{GRanges} of the plotted region.}
#'   \item{creSignalPlotGR}{\link[GenomicRanges]{GRanges} of the CRE signal
#'   (signed negative log CRE p-value) across the plotted region.}
#'   \item{assocGinter}{\link[InteractionSet]{GInteractions} of the DegCre
#'   associations in the plotted region.}
#' }
#' @examples
#' #Load example data.
#' data(DexNR3C1)
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
#'                                    DegP=DexNR3C1$DegGR$pVal,
#'                                    DegLfc=DexNR3C1$DegGR$logFC,
#'                                    CreGR=DexNR3C1$CreGR,
#'                                    CreP=DexNR3C1$CreGR$pVal,
#'                                    CreLfc=DexNR3C1$CreGR$logFC)
#'
#' #Make browser plot from specified gene name.
#' browserOuts <- plotBrowserDegCre(degCreResList=degCreResListDexNR3C1,
#'                                  geneName="ERRFI1",
#'                                  geneNameColName="GeneSymb",
#'                                  CreSignalName="NR3C1")
#' dev.off()
#'
#' #Make plot of specified region.
#' zoomGR <- GRanges(seqnames="chr1",ranges=IRanges(start=7900e3,end=8400e3))
#'
#' zoomedBrowserOuts <- plotBrowserDegCre(degCreResList=degCreResListDexNR3C1,
#'                                        plotRegionGR=zoomGR,
#'                                        geneNameColName="GeneSymb",
#'                                        CreSignalName="NR3C1")
#' dev.off()
#'
#' @author Brian S. Roberts
#'
#' @export
plotBrowserDegCre <- function(degCreResList,
                              assocAlpha = 0.05,
                              browserWinPad = 1000,
                              geneName = NULL,
                              plotRegionGR = NULL,
                              CreSignalName = "CRE",
                              assembly = "hg38",
                              plotWidth = dev.size("in")[1],
                              plotHeight = dev.size("in")[2],
                              plotXbegin = 0.9,
                              mergeGenePromotersDist = 1000,
                              sigPlotMaxY = 4,
                              assocColorRange = NULL,
                              lowAssocColor = "#88CCEE",
                              hiAssocColor = "#CC6677",
                              signalColor = "#DDCC77",
                              geneLabelFontsize = 8,
                              axisFontSize = 6,
                              panelTitleFontSize = 7,
                              geneNameColName = NULL,
                              geneHighlightDf = NULL,
                              dePrioritizeSmallRNA = TRUE,
                              useLogFC = TRUE,
                              creSignalBinRes = 100){

  if(all(c(is.null(geneName),is.null(plotRegionGR)))){
    warning("Gene name or plot region GenomicRanges must be provided")
    return(NA)
  }

  promCreHitsX <- degCreResList$degCreHits

  keepPromCreHitsX <- promCreHitsX[which(mcols(promCreHitsX)$assocProbFDR<=
    assocAlpha)]

  DegGRX <- degCreResList$DegGR
  degGrMcolsDf <- as.data.frame(mcols(DegGRX))

  CreGRX <- degCreResList$CreGR
  creGrMcolsDf <- as.data.frame(mcols(CreGRX))

  if(is.null(geneNameColName)){
    #try to find geneName in the mcols of DegGR to find the right
    #mcols columns name
    listGrepsX <- lapply(seq_len(ncol(degGrMcolsDf)),function(colX){
      geneGreppy <- paste0("^",geneName,"$")
      valColX <- degGrMcolsDf[,colX]
      return(grep(geneGreppy,valColX))
    })

    lengthsGrep <- unlist(lapply(listGrepsX,length))
    maskFoundGene <- which(lengthsGrep>0)

    #exceptions when more than entry matches geneName or no matches
    if(length(maskFoundGene)>1){
      warning("Multiple DegGR mcols columns have matches to gene name.
              Re-run setting geneNameColName to correct column name.")
      return(NA)
    }
    if(length(maskFoundGene)==0){
      warning("Cannot find gene name column in DEG GR.
              Re-run setting geneNameColName to correct column name.")
      return(NA)
    }
    else{
      geneNameColName <- colnames(degGrMcolsDf)[maskFoundGene]
    }
  }
  else{
    maskFoundGene <- which(colnames(degGrMcolsDf)==geneNameColName)
  }

  #if gene name is given, define plotRegionGR by all sig CREs associated
  #with it and by gene sig associated with them
  if(!is.null(geneName)){
    message("Using ",
      geneNameColName," as gene name column.")

    plotRegionGR <- getDegCrePlotRegionFromGene(
      degCreResList=degCreResList,
      geneName=geneName,
      geneNameColName=geneNameColName,
      assocAlpha=assocAlpha)

    if(any(is.na(plotRegionGR))){
      return(NA)
    }
  }
  else{
    geneName <- NULL
  }

  #pad plotRegionGR by browserWinPad
  plotRegionGR <- plotRegionGR + browserWinPad

  plotChrX <- as.character(seqnames(plotRegionGR))[1]
  plotStartX <- start(plotRegionGR)[1]
  plotEndX <- end(plotRegionGR)[1]

  plotUnitX <- round((plotEndX - plotStartX)/creSignalBinRes,0)

  message("Minimum display region size = ",plotUnitX," bp.")
  message("All smaller regions will be expanded.")

  #get Cre pvals in padded GR
  padCreGR <- creGRToSignal(CreGR=CreGRX,
    plotRegionGR=plotRegionGR,
    useLogFC=TRUE,
    creSignalBinRes=creSignalBinRes)

  #now make Ginteractions for associations
  #find what associations within plotRegionGR
  #associations must be in keepPromCreHitsX (pass assocAlpha)

  degIndicesInRegion <- queryHits(findOverlaps(DegGRX,plotRegionGR))

  keepDegMask <- which(queryHits(keepPromCreHitsX) %in% degIndicesInRegion)

  if(length(keepDegMask)==0){
    warning("No associations passing assocAlpha found for specified region.")
    return(NA)
  }

  GInterX <- makePlotGInter(DegGR=DegGRX,
    CreGR=CreGRX,
    plotRegionGR=plotRegionGR,
    keepPromCreHits=keepPromCreHitsX,
    plotUnit=plotUnitX,
    maskFoundGene=maskFoundGene,
    mergeGenePromotersDist=mergeGenePromotersDist)

  #make browser plot
  pageCreate(width=plotWidth,height=plotHeight, showGuides=FALSE,
    default.units = "inches",xgrid = 0, ygrid = 0)

  subPlotWidth <- plotWidth - plotXbegin - 0.1

  currentY <- 0.02

  #plot associations as arches
  archesPlotH <- 0.5
  archesAlpha <- 0.8

  assocColorPallette <- colorRampPalette(c(lowAssocColor, hiAssocColor))

  if(is.null(assocColorRange)){
    assocColorRange <- c(0,max(GInterX$assocProb))
  }

  assocColorBy <- plotgardener::colorby("assocProb",
    palette=assocColorPallette,
    range=assocColorRange)

  assocArchesPlot <- plotgardener::plotPairsArches(GInterX,
    chrom= plotChrX,
    chromstart= plotStartX,
    chromend= plotEndX,
    draw=TRUE,
    fill= assocColorBy,
    linecolor=NA,
    flip=FALSE,
    x=plotXbegin,
    y=currentY,
    width=subPlotWidth,
    height=archesPlotH,
    alpha=archesAlpha)

  #add a gradient color bar
  gradientColBarWidth <- 0.125

  colKeyX <- plotXbegin-gradientColBarWidth

  if(length(GInterX)>0){
    plotgardener::annoHeatmapLegend(plot = assocArchesPlot,
      fontcolor = "black",
      x = colKeyX,
      y = currentY+0.075*archesPlotH,
      width = gradientColBarWidth,
      height = 0.85*archesPlotH,
      fontsize = axisFontSize)
  }

  currentY <- currentY + archesPlotH +0.07

  #now make a plot of CRE differential p
  #plot signal plot
  diffSignalPlotRangeX <- c(-sigPlotMaxY,sigPlotMaxY)
  diffSignalYAnnotMax <- sigPlotMaxY

  signalPlotH <- 0.5

  crePvalplot <- plotgardener::plotSignal(data= padCreGR,
    negData=TRUE,
    chrom= plotChrX,
    chromstart= plotStartX,
    chromend= plotEndX,
    draw=TRUE,
    linecolor=c(NA,NA),
    fill=signalColor,
    x=plotXbegin,
    y=currentY,
    width=subPlotWidth,
    height=signalPlotH,
    range= diffSignalPlotRangeX,
    baseline.color="black",
    baseline.lwd=1)

  #create y axis label
  plotgardener::annoYaxis(crePvalplot,
    at=c(-diffSignalYAnnotMax,0, diffSignalYAnnotMax),
    fontsize=axisFontSize,
    axisLine=TRUE)


  currentY <- currentY + signalPlotH +0.1

  #plot genomic position
  genPosPlotH <- 0.2

  genomicPosPlot <- plotgardener::plotGenomeLabel(chrom= plotChrX,
    chromstart= plotStartX,
    chromend= plotEndX,
    x=plotXbegin,
    y=currentY,
    length=subPlotWidth,
    height=genPosPlotH,
    scale="Kb",
    assembly=assembly,
    fontsize=geneLabelFontsize)

  currentY <- currentY + genPosPlotH

  #now plot gene models
  genesPlotH <- 0.5

   testGenePlot <- plotgardener::plotGenes(chrom= plotChrX,
     chromstart= plotStartX,
     chromend= plotEndX,
     x=plotXbegin,
     y=currentY,
     width=subPlotWidth,
     height=genesPlotH,
     fontsize=geneLabelFontsize,
     draw=FALSE,
     assembly=assembly)

   testPlottedGenes <- testGenePlot$genes

   #de-prioritize plotting miRNA and snoRNA if flag is set
   if(dePrioritizeSmallRNA){
     maskSmallRNA <- grep("^MIR|^SNORD",testPlottedGenes)
     if(length(maskSmallRNA)>0){
       smallRNAGeneNames <- testPlottedGenes[maskSmallRNA]
       testPlottedGenes <-
        c(setdiff(testPlottedGenes,smallRNAGeneNames),smallRNAGeneNames)
     }
   }

   if(!is.null(geneName) & is.null(geneHighlightDf)){
     genePlotOrder <- c(geneName,setdiff(testPlottedGenes,geneName))
     geneHighlightDf <- data.frame(gene=genePlotOrder,
       color=c("black",rep("darkgray",
       (length(setdiff(testPlottedGenes,geneName))))))
   }
   else{
     if(!is.null(geneHighlightDf)){
       genePlotOrder <- testPlottedGenes
       geneHighlightDf <- geneHighlightDf
     }
     else{
       genePlotOrder <- testPlottedGenes
       geneHighlightDf <- data.frame(gene=genePlotOrder,color="black")
     }
   }

   genePlot <- plotgardener::plotGenes(chrom= plotChrX,
     chromstart= plotStartX,
     chromend= plotEndX,
     x=plotXbegin,
     y=currentY,
     width=subPlotWidth,
     height=genesPlotH,
     fontsize=geneLabelFontsize,
     draw=TRUE,
     assembly=assembly,
     geneOrder=genePlotOrder,
     geneHighlights=geneHighlightDf,)

   #add panel titles
   panelTitleX <- 0.4*plotXbegin

   #arches plot
   archesPlotY <- getLabelYfromPlotgardenerObj(assocArchesPlot)
   archesPlotLabel <- "DegCre\nAssociation\nProbability"

   plotgardener::plotText(label=archesPlotLabel,
    fontsize=panelTitleFontSize,
    x=panelTitleX,
    y=archesPlotY,
    just="center")

  #signal plot
  creSignaPlotY <- getLabelYfromPlotgardenerObj(crePvalplot)

   if(useLogFC){
     creSignalPlotLabel1 <- paste(CreSignalName,"signed")
   }
   else{
     creSignalPlotLabel1 <- CreSignalName
   }

   creSignalPlotLabel2 <- bquote("-log"[10]~"(p)")

   plotgardener::plotText(label=creSignalPlotLabel1,
    fontsize=panelTitleFontSize,
    x=panelTitleX,
    y=as.numeric(creSignaPlotY)-0.05*signalPlotH,
    just="bottom")

  plotgardener::plotText(label=creSignalPlotLabel2,
    fontsize=panelTitleFontSize,
    x=panelTitleX,
    y=as.numeric(creSignaPlotY)+0.05*signalPlotH,
    just="top")

  outList <- list()
  outList$plotRegionGR <- plotRegionGR
  outList$creSignalPlotGR <- padCreGR
  outList$assocGinter <- GInterX

  invisible(outList)
}


#' Convert CreGR to Pseudo-Continuous Signal for Plotting
#'
#' This function converts a GenomicRanges object containing CRE data
#' to a pseudo-continuous signal track suitable for plotting in plotgardener.
#' The signal is derived from p-values and, optionally, log-fold change values
#' associated with CREs.
#'
#' @param CreGR A \link[GenomicRanges]{GRanges} object representing CRE data.
#' @param plotRegionGR A \link[GenomicRanges]{GRanges} object specifying the
#' region of interest for plotting.
#' @param useLogFC Logical, indicating whether to consider log-fold change
#' values (Default: \code{TRUE}).
#' @param pValColName Character specifying the column name in CreGR containing
#' p-values (Default: \code{pVal}).
#' @param logFcColName Character specifying the column name in CreGR
#' containing log-fold change values (Default: \code{logFC}).
#' @param creSignalBinRes Numeric value specifying the bin resolution in base
#' pairs for the pseudo-continuous signal (Default: \code{100}).
#'
#' @return A \link[GenomicRanges]{GRanges} object with signal values in
#' metadata column \code{score}, suitable for plotting in
#' \link[plotgardener]{plotSignal}.
#'
#' @details
#' Not exported. This function takes a \link[GenomicRanges]{GRanges} object
#' (\code{CreGR}) representing CRE data, extracts p-values, and,
#' if specified, log-fold change values. It then converts these values into a
#' signed -log p-value in the \code{signal} metadata column.
#'
#' @keywords internal
#'
#' @examples
#' #Load example data.
#' data(DexNR3C1)
#'
#' myCreGR <- DexNR3C1$CreGR
#' myPlotRegionGR <- GRanges(seqnames="chr4",
#'                           ranges=IRanges(start=3.6e6,end=3.8e6))
#'
#' # Convert CRE data to a pseudo-continuous signal
#' creSignalGR <- creGRToSignal(CreGR=myCreGR, plotRegionGR=myPlotRegionGR)
#'
#' @author Brian S. Roberts
#'
creGRToSignal <- function(
  CreGR,
  plotRegionGR,
  useLogFC=TRUE,
  pValColName="pVal",
  logFcColName="logFC",
  creSignalBinRes=100){

  creGrMcolsDf <- mcols(CreGR)
  rawCrePVals <- creGrMcolsDf[[pValColName]]

  #get logFCs
  if(!useLogFC){
    logFcs <- rep(1,nrow(creGrMcolsDf))
  }
  else{
    logFcs <- creGrMcolsDf[[logFcColName]]
  }

  signedCreNegLog10PVals <- -1*log10(rawCrePVals)*sign(logFcs)
  maskNoNA <- which(!is.na(signedCreNegLog10PVals))
  noNaNegLog10PVals <- signedCreNegLog10PVals[maskNoNA]

  #pad Cre signal for nice plotting
  padCreGR <- unlist(tile(plotRegionGR+10,n=creSignalBinRes))
  padCreGR$score <- rnorm(n=creSignalBinRes,mean=0,sd=0.0001)

  noNaCreGR <- CreGR[maskNoNA]

  hitsPadtoCre <- findOverlaps(padCreGR,noNaCreGR)

  realCreScores <- tapply(noNaNegLog10PVals[subjectHits(hitsPadtoCre)],
    INDEX=queryHits(hitsPadtoCre),FUN=function(x){

    maskHighestAbs <- which(abs(x)==max(abs(x)))
    if(length(maskHighestAbs)>1){
      maskHighestAbs <- maskHighestAbs[1]
    }
    outVal <- x[maskHighestAbs]
    return(outVal)
  })

  padCreGR$score[unique(queryHits(hitsPadtoCre))] <- realCreScores
  return(padCreGR)
}


#' Create a GInteractions Object for Plotgardener Arches Plotting
#'
#' This function generates a GInteractions object suitable for plotting arches
#' in plotgardener.
#'
#' @param DegGR A \link[GenomicRanges]{GRanges} object representing DEG TSSs.
#' @param CreGR A \link[GenomicRanges]{GRanges} object representing CRE regions.
#' @param plotRegionGR A \link[GenomicRanges]{GRanges} object specifying the
#' region of interest for plotting.
#' @param keepPromCreHits A \link[S4Vectors]{Hits} object containing filtered
#' promoter-CRE associations.
#' @param plotUnit Numeric value specifying the minimum plot unit width in
#' base pairs.
#' @param maskFoundGene Character specifying the index of column for gene name
#' in \code{DegGR}.
#' @param mergeGenePromotersDist Numeric value specifying the distance in base
#' pairs within which promoters of the same gene are merged.
#'
#' @return A \link[InteractionSet]{GInteractions} object suitable for arches
#' plotting in plotgardener.
#'
#' @details
#' This function takes genomic data, including differential expression data
#' (\code{DegGR}), CRE data (\code{CreGR}),
#' a region of interest (\code{plotRegionGR}), TSS-CRE associations
#' (\code{keepPromCreHits}), and other parameters,
#' and creates a GInteractions object suitable for plotting arches in
#' plotgardener.
#'
#' Not exported. The function processes the data to ensure that promoters and
#' CREs are in a suitable format for plotting. It merges
#' promoters within a specified distance (\code{mergeGenePromotersDist}) and
#' adjusts the width of regions to meet the
#' minimum plot unit width (\code{plotUnit}). It also handles cases where
#' promoters and CREs may have been merged or
#' contain duplicate regions, ensuring correct associations.
#' This function is meant to run within \link{plotBrowserDegCre}. It will not
#' run well on unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Create the GInteractions object for plotting
#' GInterX <- makePlotGInter(DegGR=myDegGR,
#'                           CreGR=myCreGR,
#'                           plotRegionGR=myPlotRegionGR,
#'                           keepPromCreHits=myKeepPromCreHits,
#'                           plotUnit=5000,
#'                           maskFoundGene=3,
#'                           mergeGenePromotersDist=2000)
#' }
#'
#' @author Brian S. Roberts
#'
makePlotGInter <- function(DegGR,
                         CreGR,
                         plotRegionGR,
                         keepPromCreHits,
                         plotUnit,
                         maskFoundGene,
                         mergeGenePromotersDist){


  degIndicesInRegion <- queryHits(findOverlaps(DegGR,plotRegionGR))
  creIndicesInRegion <- queryHits(findOverlaps(CreGR,plotRegionGR))

  keepDegMask <- which(queryHits(keepPromCreHits) %in% degIndicesInRegion)
  keepCreMask <- which(subjectHits(keepPromCreHits) %in% creIndicesInRegion)

  keepAllMask <- sort(unique(c(keepDegMask,keepCreMask)))

  inRegionPromCreHits <- keepPromCreHits[keepAllMask]

  #merge promoters for same gene that are within mergeGenePromotersDist
  rawInRegionDegGR <- DegGR[queryHits(inRegionPromCreHits)]
  strand(rawInRegionDegGR) <- rep("*",length(rawInRegionDegGR))

  rawInRegionDegMcolsDf <- as.data.frame(mcols(rawInRegionDegGR))
  inRegionGeneNames <- as.character(rawInRegionDegMcolsDf[,maskFoundGene])
  uniqInRegionGeneNames <- unique(inRegionGeneNames)
  names(uniqInRegionGeneNames) <- uniqInRegionGeneNames

  listRawInRegionMergedDegs <- lapply(uniqInRegionGeneNames,
    function(uniqGeneX){
      maskX <- which(inRegionGeneNames == uniqGeneX)
      inRegionDegGR <- rawInRegionDegGR[maskX]
      if(length(maskX) == 1){
        outGRX <- granges(inRegionDegGR)
        outGRX$origI <- maskX
      }
      else{
        tempMergeGRx <-
          GenomicRanges::reduce(inRegionDegGR+mergeGenePromotersDist)
        hitsXtoTempMerge <-
          GenomicRanges::findOverlaps(tempMergeGRx,inRegionDegGR)
        DegGRsubjHits <- inRegionDegGR[subjectHits(hitsXtoTempMerge)]
        outGRx <- tapply(GenomicRanges::granges(DegGRsubjHits),
          INDEX=queryHits(hitsXtoTempMerge),FUN=function(setGRx){
            return(range(setGRx))
        })
        mergeGRx <- do.call(c,unname(outGRx))
        if(length(mergeGRx) < length(inRegionDegGR)){
          message(uniqGeneX,
            " promoters have been merged for plotting.")
        }
        hitsInRegiontoMerged <-
          GenomicRanges::findOverlaps(inRegionDegGR,mergeGRx)
        outGRX <- mergeGRx[subjectHits(hitsInRegiontoMerged)]
        outGRX$origI <- maskX
      }

      return(outGRX)
    })

  rawInRegionMergedDegGR <- do.call(c,unname(listRawInRegionMergedDegs))

  tempInRegionDegGR <- rawInRegionMergedDegGR[
    order(rawInRegionMergedDegGR$origI)]

  tempInRegionDegGR <- unique(tempInRegionDegGR)

  #make the DegGR regions a minimum width for pretty plotting
  maskDegTooSmall <- which(width(tempInRegionDegGR)<plotUnit)

  if(length(maskDegTooSmall)>0){
    tempInRegionDegGR[maskDegTooSmall] <-
      tempInRegionDegGR[maskDegTooSmall]*
      (width(tempInRegionDegGR[maskDegTooSmall])/plotUnit)
  }


  #make the CreGR regions a minimum width for pretty plotting
  rawInRegionCreGR <- CreGR[subjectHits(inRegionPromCreHits)]

  tempInRegionCreGR <- granges(rawInRegionCreGR)

  strand(tempInRegionCreGR) <- rep("*",length(tempInRegionCreGR))

  #unique regions, may have repeats
  tempInRegionCreGR <- unique(tempInRegionCreGR)

  maskCreTooSmall <- which(width(tempInRegionCreGR)<plotUnit)

  if(length(maskCreTooSmall)>0){
    tempInRegionCreGR[maskCreTooSmall] <-
      tempInRegionCreGR[maskCreTooSmall]*
      (width(tempInRegionCreGR[maskCreTooSmall])/plotUnit)
  }

  #expanding tempInRegionCreGR for plot readability could yield to
  #overlapping regions which will cause problems later. Reduce result
  #to avoid this
  tempInRegionCreGR <- reduce(tempInRegionCreGR)

  #both tempInRegionDegGR and tempInRegionCreGR may have been merged
  #or uniqued. Will no longer match up in inRegionPromCreHits. Need
  #to map back to original to get correct association structure

  #deg
  hitsTempToRawInRegionDegGR <-
    GenomicRanges::findOverlaps(tempInRegionDegGR,
    rawInRegionDegGR,maxgap=mergeGenePromotersDist)

  origRawDegIndices <-
    queryHits(inRegionPromCreHits)[subjectHits(hitsTempToRawInRegionDegGR)]

  maskDegTempToRaw <- match(queryHits(inRegionPromCreHits),origRawDegIndices)

  newTempDegQueryHits <-
    queryHits(hitsTempToRawInRegionDegGR)[maskDegTempToRaw]

  #cre
  hitsTempToRawInRegionCreGR <- GenomicRanges::findOverlaps(tempInRegionCreGR,
    rawInRegionCreGR)

  origRawCreIndices <-
    subjectHits(inRegionPromCreHits)[
      subjectHits(hitsTempToRawInRegionCreGR)]

  maskCreTempToRaw <-
    match(subjectHits(inRegionPromCreHits),origRawCreIndices)

  newTempCreSubjHits <-
    queryHits(hitsTempToRawInRegionCreGR)[maskCreTempToRaw]

  newTempRefHitsMat <- cbind(newTempDegQueryHits,newTempCreSubjHits)

  #newTempRefHitsMat is analogous to inRegionPromCreHits, except that
  #the indices refer to tempInRegionDegGR and tempInRegionCreGR
  #rather than DegGR and CreGR

  #newTempRefHitsMat may contain replicate rows that will now be squeezed
  hashNewTempRefHitsMat <- paste(newTempDegQueryHits,newTempCreSubjHits,
    sep="_")

  uniqHashNewTempRefHitsMat <- unique(hashNewTempRefHitsMat)

  if(length(uniqHashNewTempRefHitsMat)<length(hashNewTempRefHitsMat)){
    listNewTempRefHitsMat <- tapply(seq_len(nrow(newTempRefHitsMat)),
      INDEX=hashNewTempRefHitsMat,FUN=function(iX){
        return(newTempRefHitsMat[iX[1],])
      })
    newTempRefHitsMat <- matrix(unlist(listNewTempRefHitsMat),
      byrow=TRUE,ncol=2)

    newAdjAssocProb <- tapply(mcols(inRegionPromCreHits)$assocProb,
      INDEX=hashNewTempRefHitsMat,function(x){
        return(max(x,na.rm=TRUE))
      })

    newAdjAssocProb <- as.numeric(newAdjAssocProb)

    newAdjAssocProbFDR <- tapply(mcols(inRegionPromCreHits)$assocProbFDR,
      INDEX=hashNewTempRefHitsMat,function(x){
        return(min(x,na.rm=TRUE))
      })

    newAdjAssocProbFDR <- as.numeric(newAdjAssocProbFDR)
  }
  else{
    newAdjAssocProb <- mcols(inRegionPromCreHits)$assocProb
    newAdjAssocProbFDR <- mcols(inRegionPromCreHits)$assocProbFDR
  }

  subTempInRegionDegGR <-
    GenomicRanges::granges(tempInRegionDegGR[newTempRefHitsMat[,1]])

  subTempInRegionCreGR <-
    GenomicRanges::granges(tempInRegionCreGR[newTempRefHitsMat[,2]])

  GInterX <-
    InteractionSet::GInteractions(subTempInRegionDegGR,
      subTempInRegionCreGR)

  GInterX$assocProb <- newAdjAssocProb
  GInterX$assocProbFDR <- newAdjAssocProbFDR
  return(GInterX)
}


#' Get Genomic Range for a Gene and Associated CREs Below an FDR Threshold
#'
#' Given a list of DegCre results, (\code{degCreResList}), this function
#' generates a \link[GenomicRanges]{GRanges} object
#' encompassing all associated CRE regions for a specific gene with an
#' associated FDR below a specified threshold.
#'
#' @param degCreResList List of DegCre results.
#' @param geneName Character of the name of the gene for which to retrieve
#' associated CRE regions.
#' @param geneNameColName Character specifying the column name for gene names
#' in \code{DegGR}.
#' @param assocAlpha Numeric value from 0 to 1 specifying the threshold for the
#' association probability FDR (Default: \code{0.05}).
#'
#' @return A \link[GenomicRanges]{GRanges} object representing the genomic
#' region encompassing all associated CREs for the specified gene,
#'   or \code{NA} if no associations below the FDR threshold are found.
#'
#' @details
#' Not exported. This function extracts the relevant components from the input
#' \code{degCreResList} and identifies associations for the
#' specified gene with an  association probability FDR below \code{assocAlpha}.
#' If associations are found, it computes
#' the genomic range encompassing all associated CREs and returns it as a
#' GenomicRanges object. If no associations meet the
#' threshold, it returns \code{NA}.
#' This function is meant to run within \link{plotBrowserDegCre}. It will not
#' run well on unintended inputs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' #Load example data.
#' data(DexNR3C1)
#'
#' #Generate DegCre results.
#' degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
#'                          DegP=DexNR3C1$DegGR$pVal,
#'                          DegLfc=DexNR3C1$DegGR$logFC,
#'                          CreGR=DexNR3C1$CreGR,
#'                          CreP=DexNR3C1$CreGR$pVal,
#'                          CreLfc=DexNR3C1$CreGR$logFC)
#'
#' #Get plot region.
#' plotRegionGR <-
#'  getDegCrePlotRegionFromGene(degCreResList=degCreResListDexNR3C1,
#'                              geneName = "ERRFI1",
#'                              geneNameColName = "GeneSymb",
#'                              assocAlpha = 0.05)
#' }
#'
#' @author Brian S. Roberts
#'
getDegCrePlotRegionFromGene <- function(degCreResList,
                                      geneName,
                                      geneNameColName,
                                      assocAlpha=0.05){

  DegGRX <- degCreResList$DegGR
  CreGRX <- degCreResList$CreGR

  degCreHits <- degCreResList$degCreHits
  keepDegCreHits <- degCreHits[which(mcols(degCreHits)$assocProbFDR<=
    assocAlpha)]

  degGrMcolsDf <- as.data.frame(mcols(DegGRX))

  geneNameCol <- which(colnames(degGrMcolsDf)==geneNameColName)

  maskGeneX <- which(degGrMcolsDf[,geneNameCol] == geneName)
  maskGeneInHits <- which(queryHits(keepDegCreHits) %in% maskGeneX)

  if(length(maskGeneInHits)==0){
    warning("No associations passing assocAlpha found for ",geneName)
    plotRegionGR <- NA
  }

  allFoundSubJHits <- subjectHits(keepDegCreHits[maskGeneInHits])
  maskHitsSubj <- which(subjectHits(keepDegCreHits) %in%
    allFoundSubJHits)

  geneXPromCreHits <- keepDegCreHits[maskHitsSubj]

  allDegCreGRX <-
    c(GenomicRanges::granges(DegGRX[queryHits(geneXPromCreHits)]),
    GenomicRanges::granges(CreGRX[subjectHits(geneXPromCreHits)]))
  plotRegionGR <- range(allDegCreGRX,ignore.strand=TRUE)

  return(plotRegionGR)
}


#' Get Y Coordinate for Label Placement from a PlotGardener Plot Object
#'
#' Given a PlotGardener plot object, this function calculates the Y coordinate
#' for label placement based on the vertical positioning specified by the
#' 'just' parameter.
#'
#' @param plotgardenerObj A plotgardener plot object.
#' @param just Character specifying the vertical positioning of the label.
#'   Options include "top" (top-aligned), "bottom" (bottom-aligned), and
#'   "center" (center-aligned). (Default: \code{center})
#'
#' @return A numeric value representing the Y coordinate in inches for label
#'  placement.
#'
#' @details
#' Not exported. This function takes a plotgardener plot object as input and
#' calculates the Y coordinate for label placement based on the vertical
#' positioning specified by the 'just' parameter.
#' The 'just' parameter determines whether the label should be placed at the
#' top, bottom, or center of the plot object. The function returns the
#' calculated Y coordinate.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get y label.
#' labelY <- getLabelYfromPlotgardenerObj(plotObj=mySignalPlotObj,
#'                                        just = "center")
#' }
#'
#' @author [Author Name]
#'
getLabelYfromPlotgardenerObj <- function(plotgardenerObj,
                                       just="center"){

  plotY <- as.numeric(plotgardenerObj$y)
  plotH <- as.numeric(plotgardenerObj$height)
  if(just=="top"){
    outY <- plotY
  }
  if(just=="bottom"){
    outY <- plotY + plotH
  }
  if(just=="center"){
    outY <- plotY + 0.5*plotH
  }
  return(outY)
}
