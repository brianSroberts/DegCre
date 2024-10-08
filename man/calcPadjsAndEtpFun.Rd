% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degCreCore.R
\name{calcPadjsAndEtpFun}
\alias{calcPadjsAndEtpFun}
\title{Convert p-values to adjusted p-values}
\usage{
calcPadjsAndEtpFun(pVal, method = "qvalue")
}
\arguments{
\item{pVal}{Input p-values.}

\item{method}{One of "bonferroni" or "qvalue". (Default: \code{qvalue}).}
}
\value{
A list two slots:
\describe{
\item{pAdj}{Adjusted p-values.}
\item{etpFun}{A function for calculating the number of expected true
 positives from adjusted p-values.}
}
}
\description{
Convert a set of p-values to adjusted pvalues with a specified method and
return a function for converting set of pvalues to number of 
expected true positives.
}
\details{
Not exported. This function adjusted p-values by the specified method. It
also returns a function that calculates the fraction of expected true
positives from a subset of the adjusted p-values
}
\examples{
\dontrun{
outPadjList <- calcPadjsAndEtpFun(myPVals,method="qvalue")
}
}
\author{
Brian S. Roberts
}
\keyword{internal}
