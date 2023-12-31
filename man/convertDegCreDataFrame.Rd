% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degCreConversions.R
\name{convertDegCreDataFrame}
\alias{convertDegCreDataFrame}
\title{Convert DegCre Results List to DataFrame}
\usage{
convertDegCreDataFrame(degCreResList, assocAlpha = 0.05)
}
\arguments{
\item{degCreResList}{List of DegCre results.}

\item{assocAlpha}{The significance threshold for associations to be included
in the output (Default: \code{0.05}).}
}
\value{
A \link[S4Vectors]{DataFrame} containing the significant associations
that pass the specified significance threshold. It is roughly in BEDPE
format.
}
\description{
Given a DegCre results list, this function converts it into a DataFrame for
further analysis and export.
}
\details{
This function takes a DegCre results list as input and extracts the
significant associations based on the adjusted p-values \code{assocProbFDR}
compared to the specified significance threshold \code{assocAlpha}.
It then creates a \link[S4Vectors]{DataFrame} with the genomic coordinates
of the significant associations from both the \code{DegGR} and \code{CreGR}
components of the input list.
These are marked as \code{Deg_} or \code{Cre_} with \code{chr}, \code{start},
\code{end}, and \code{strand}.
The coordinates are followed by the metadata of the \link[S4Vectors]{Hits}
\link[S4Vectors]{DataFrame}d by \link{runDegCre}.
These are then followed by all metadata columns in the input \code{DegGR} or
\code{CreGR} proceeded by either \code{Deg_} or \code{Cre_} in the colname.

If no associations pass the significance threshold,
the function returns \code{NA}.
}
\examples{
#Load test data.
data(DexNR3C1)

#Generate DegCre results.
degCreResListDexNR3C1 <- runDegCre(DegGR=DexNR3C1$DegGR,
                                   DegP=DexNR3C1$DegGR$pVal,
                                   DegLfc=DexNR3C1$DegGR$logFC,
                                   CreGR=DexNR3C1$CreGR,
                                   CreP=DexNR3C1$CreGR$pVal,
                                   CreLfc=DexNR3C1$CreGR$logFC)

#Create DataFrame.
outDf <-
  convertDegCreDataFrame(degCreResList=degCreResListDexNR3C1,
                                    assocAlpha = 0.05)

#Write out as text file.

write.table(outDf,file="myDegCreResults.tsv",sep="\t",row.names=F,quote=F)

}
\author{
Brian S. Roberts
}
