% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapseDegCreToGene.R
\name{collapseDegCreToGene}
\alias{collapseDegCreToGene}
\title{Collapse DegCre associations to a single TSS per gene}
\usage{
collapseDegCreToGene(
  degCreResList,
  method = "nearest",
  geneColname = "GeneSymb",
  useParallel = TRUE
)
}
\arguments{
\item{degCreResList}{List of DegCre results.}

\item{method}{Method for choosing between multiple TSS. Currently only
supported for "nearest". (Default:\code{"nearest"})}

\item{geneColname}{The name of the metadata column in DegGR within
\code{degCreResList} that has the gene name (Default:\code{"GeneSymb"})}

\item{useParallel}{Logical whether to implement parallelization with
\link[BiocParallel]{BiocParallel} package. (Default:\code{FALSE})}
}
\value{
A \code{degCreResList} with only the shortest TSS to CRE
association per gene.
}
\description{
Given a DegCre results list, this function finds associations between
the same CRE and multiple TSSs of the same gene and keeps the nearest
TSS only.
}
\details{
Often, the DegGR input to DegCre will contain multiple TSS's for a given
gene. DegCre will create associations to all of them. Often, downstream
analyses need to only have on CRE to TSS association per gene. This function
modifies the DegCre \link[S4Vectors]{Hits} to keep the shortest (minimum)
genomic distance association.
}
\examples{
#Load required packages.
library(GenomicRanges)

#Load test data.
data(DexNR3C1)

subDegGR <-
 DexNR3C1$DegGR[which(GenomeInfoDb::seqnames(DexNR3C1$DegGR)=="chr1")]
subCreGR <-
 DexNR3C1$CreGR[which(GenomeInfoDb::seqnames(DexNR3C1$CreGR)=="chr1")]

#Generate DegCre results.
degCreResListDexNR3C1 <- runDegCre(DegGR=subDegGR,
                                   DegP=subDegGR$pVal,
                                   DegLfc=subDegGR$logFC,
                                   CreGR=subCreGR,
                                   CreP=subCreGR$pVal,
                                   CreLfc=subCreGR$logFC)

#Convert to single TSS per association.

degCreResListUniqTSS <- collapseDegCreToGene(degCreResListDexNR3C1,
                                             method = "nearest",
                                             geneColname = "GeneSymb",
                                             useParallel=FALSE)

}
\author{
Brian S. Roberts
}
