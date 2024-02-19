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
