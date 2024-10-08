% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degCreCore.R
\name{makeDistBinChunkList}
\alias{makeDistBinChunkList}
\title{Convert a Hits object to association distance list}
\usage{
makeDistBinChunkList(sortHitsWithDistDf, pickedBinSize, verbose = FALSE)
}
\arguments{
\item{sortHitsWithDistDf}{A \link[S4Vectors]{DataFrame} with queryHits and
subjectHits from a \link[S4Vectors]{Hits} object. They must be sorted by
increasing association distance.}

\item{pickedBinSize}{The number of associations per bin.}
}
\value{
A list of \link[S4Vectors]{DataFrame}(s) with 
\code{pickedBinSize} associations, grouped by distance.
}
\description{
Convert a Hits of DEG to CRES with distance to a list broken into distance
bins
}
\details{
Not exported. This function splits by the associations by distance to create
bins with equal number of elements. The last bin is merged to the second to
last bin if it is smaller than 0.8*\code{pickedBinSize}.
}
\examples{
\dontrun{
outChunkList <- makeDistBinChunkList(inSortHitsWithDistDf,pickedBinSize=1000)
}
}
\author{
Brian S. Roberts
}
\keyword{internal}
