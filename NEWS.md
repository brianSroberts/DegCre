# CHANGES IN VERSION 1.2.0
* Version up to stable release.

# CHANGES IN VERSION 1.1.0.9000

## NEW FEATURES

* `runDegCre` now accepts `"qvalue"` as an input to `pAdjMethod` which is now the default. Our testing shows that this method improves the performance of DegCre predictions. The qvalue calculation does not require the user to specify a `alphaVal`.

* Added new function `collapseDegCreToGene` which converts DegCre results for gene with multiple TSSs to use only the association that spans the shortest distance.

* Added new function `calcAssocProbOR` that calculates the odds-ratio for DegCre association probabilities. This function can operate on `assocProb` or `rawAssocProb` values by altering the `type` flag. It is meant to replace `calcRawAssocProbOR`.

* Added new function `convDegCreResListToCreGeneScoreGR` that converts DegCre results to a simplified `GRanges` with the predicted gene and score as metadata.

## OTHER CHANGES

* Streamlined `runDegCre` code to use more subfunctions.

## BUG FIXES

* Changed the internal function `makePlotGInter` to return a sorted `GInteractions` to avoid inconsistent returns.

# CHANGES IN VERSION 1.0.2

* Changed citation to Genome Research article

## BUG FIXES

* Added parameter `minNDegs` to `optimizeAlphaDegCre` to stop it from trying `testedAlphaVals` that result in the number of passing DEGs to be too low for the optimization algorithm to function properly.

# CHANGES IN VERSION 0.99.16

* Added utils to DESCRIPTION

# CHANGES IN VERSION 0.99.15

* Updated to use R 4.4.0

* Added importFrom for grDevices and utils

* Fixed typos in help files

# CHANGES IN VERSION 0.99.14
## BUG FIXES

* Added imported packages to DESCRIPTION

# CHANGES IN VERSION 0.99.13

* Numerous code style changes made to comply with Bioconductor standards.

* Updated to run require R 4.3.3

* Vignette style changed to Bioconductor

* Selective imports implemented

# CHANGES IN VERSION 0.99.12
## BUG FIXES

* Fixed one more bad href in help files.

# CHANGES IN VERSION 0.99.11
## BUG FIXES

* Fixed bad hrefs in help files.

# CHANGES IN VERSION 0.99.10
## BUG FIXES

* Fixed example for runDegCre.

# CHANGES IN VERSION 0.99.9
## BUG FIXES

* Fixed test functions to subset GRanges.

# CHANGES IN VERSION 0.99.8
## BUG FIXES

* Fixed runnable examples to load GenomicRanges.

# CHANGES IN VERSION 0.99.7
## BUG FIXES

* Changed test functions to run on chr1 only for efficiency.

# CHANGES IN VERSION 0.99.6
## BUG FIXES

* Changed examples to run on chr1 for computational efficiency.

# CHANGES IN VERSION 0.99.4
## BUG FIXES

* Fixed F to FALSE in example

# CHANGES IN VERSION 0.99.3
## BUG FIXES

* Set all non-exported functions to not run examples.
* Fixed vignette code.

## NEW FEATURES

* `getAssocDistHits` is now exported.

# CHANGES IN VERSION 0.99.2
## BUG FIXES

* Fixed calcAUC.Rd to not run example.

# CHANGES IN VERSION 0.99.1
## BUG FIXES

* Fixed `calcBinomFDRperBin` to give FDR = 1 for 0 assocProbs in top 10 CRE p-values.

# CHANGES IN VERSION 0.99.0
## BUG FIXES

* Added missing unit test for `optimizeAlphaDegCre`

# CHANGES IN VERSION 0.3.0
## NEW FEATURES

* Added unit tests for most functions.

## BUG FIXES

* Fixed row names in `plotDegCreAssocProbVsDist`
* Changed left margin handling and arched plot color bar placement in 
`plotBrowserDegCre`
* Added internal keyword to un-exported functions.

# CHANGES IN VERSION 0.2.0
## NEW FEATURES

* Added a `NEWS` file to track changes.
* Made numerous style changes to code in preparation for Bioconductor submission.
