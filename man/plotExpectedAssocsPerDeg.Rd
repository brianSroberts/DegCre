% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degCrePerDegFuns.R
\name{plotExpectedAssocsPerDeg}
\alias{plotExpectedAssocsPerDeg}
\title{Plot Histogram of Expected Associations per DEG}
\usage{
plotExpectedAssocsPerDeg(
  expectAssocPerDegDf,
  barOutlineColor = "#88CCEE",
  barFillColor = NULL,
  extraText = FALSE
)
}
\arguments{
\item{expectAssocPerDegDf}{\link[S4Vectors]{DataFrame} output of
\link{getExpectAssocPerDEG}.}

\item{barOutlineColor}{Color for the outline of the histogram bars.
(Default: \code{#88CCEE})}

\item{barFillColor}{Fill color for the histogram bars. If \code{NULL},
it will be derived from \code{barOutlineColor} with adjusted transparency.}

\item{extraText}{Logical, indicating whether additional text information
(Details) should be added to the plot.}
}
\value{
Invisibly, the median expected associations per DEG.
}
\description{
Plots a histogram of the expected number of associations per DEG
(Differentially Expressed Gene) based on DegCre analysis.
}
\details{
This function generates a histogram of the expected number of associations
per DEG and optionally adds additional text information to the plot, such
as DEG FDR, association FDR, and the fraction of DEGs with at least one
association.
Plot displays a dashed line a value indicating the median expected
associations per DEG.
}
\examples{
#Load example data.
data(DexNR3C1)

#Generate DegCre results.
degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
                                   DegP=DexNR3C1$DegGR$pVal,
                                   DegLfc=DexNR3C1$DegGR$logFC,
                                   CreGR=DexNR3C1$CreGR,
                                   CreP=DexNR3C1$CreGR$pVal,
                                   CreLfc=DexNR3C1$CreGR$logFC)

# Generate data frame of expected associations per DEG
expectAssocPerDegDf <-
 getExpectAssocPerDEG(degCreResList = degCreResListDexNR3C1,
                      geneNameColName = "GeneSymb",
                      assocAlpha = 0.05)

# Plot histogram of expected associations per DEG
medianExpAssocs <- plotExpectedAssocsPerDeg(expectAssocPerDegDf,
                                            barOutlineColor = "blue",
                                            extraText = TRUE)

}
\author{
Brian S. Roberts
}
