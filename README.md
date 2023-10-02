# DegCre
DegCre is an R package for the association of DEGs to CREs

DegCre associates differentially expressed genes (DEGs) with cis-regulatory elements (CREs) in a probabilistic, non-parametric approach. DegCre is intended to be applied on differential expression and regulatory signal data derived from responses to perturbations such as drug or natural ligand treatmnents. Examples include activation of nuclear receptors with pharmacologic modulators or induction of macorphages with cyotkines.

DegCre uses the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) framework for handling genomic regions and some calculations. As one input, `DegGR`, users generate differntial expression statistics for genes with methods such as [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). These
values should then be associated with gene TSSs such as those available from [EPDNew](https://epd.expasy.org/epd/) in a `GRanges`. The second input, `CreGR`, is differential regulatory signal data (in the example, NR3C1 ChIP-seq data) associated with genomic regions in a `GRanges` such as those generated by [csaw](https://bioconductor.org/packages/release/bioc/html/csaw.html). 

A complete description of the mathematical basis of the DegCre core algorithms is provided in (CITE Roberts et al). DegCre generates a `Hits` object of all associations between `DegGR` and `CreGR` within a specified maximum distance.
Associations are then binned by TSS-to-CRE distance according to an algorithm that balances resolution (many bins with few members)
versus minimization of the deviance of each bin's CRE p-value distribution from the global distribution, selecting an optimal bin size.

Next, DegCre applies a non-parametric algorithm to find concordance between DEG and CRE differential effects within bins and derives an association probability.
For all association probabilities involving one given CRE, the probabilities are adjusted to favor associations across shorter distances.
An FDR of the association probability is then estimated. Results are returned in list containing a `Hits` object and both input `GRanges`.
