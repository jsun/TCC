\name{TCC}
\docType{package}
\alias{TCC}
\alias{TCC-package}
\title{A package for differential expression analysis from tag
 count data with robust normalization strategies}
\description{
This package performs differential expression analysis from transcriptome 
data that are produced from high-throughput sequencing (HTS) and microarray
technologies. A notable feature of this 
package is to provide robust normalization methods whose strategy 
is to remove data assigned as potential differentially expressed 
genes (DEGs) before performing normalization for RNA-seq count data
(Kadota et al., 2012; Sun et al., 2013).
}
\details{
TCC is a package for differential expression analysis from 
transcriptome data produced from RNA-seq and microarray data.
This package implements some functions for calculating normalization factors,
identifying DEGs, depicting so-called M-A plot, and generating simulation data.

To utilize this package, the count matrix coupled with label information
should be stored to a \link{TCC-class} object using the \code{new} method.
All functions, except for some legacy functions
(i.e., \code{\link{NBsample}}, \code{\link{do_TbT}},
\code{\link{exactTestafterTbT}}, and \code{\link{MAplot}} for RNA-seq data)
and for two recently added functions (i.e., \code{\link{ROKU}} and
\code{\link{WAD}}) for microarray data,
used in this package require this \link{TCC-class} object.
Using this object, the \code{\link{calcNormFactors}} function calculates
normalization factors and the \code{\link{estimateDE}} function estimates
the degree of differential expression (DE) for individual genes. 
The estimated normalization factors obtained by using the
\code{\link{calcNormFactors}} function are used within the statistical
model for differential analysis in the \code{\link{estimateDE}} function.
Both two functions internally call functions from other packages
(edgeR, DESeq, baySeq, EBSeq, and samr) when specified.
TCC also provides some useful functions: \code{\link{simulateReadCounts}}
for generating simulation data with various experimental designs,
\code{\link{plot}} for depicting a M-A plot,
\code{\link{plotFCPseudocolor}} for depicting a pseudo-color image of
simulation condition that the user specified,
\code{\link{WAD}} for identifying DEGs from two-group microarray data
(single-factor design), and \code{\link{ROKU}} for identifying
tissue-specific genes from microarray data for many tissues.
}
\examples{
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
show(tcc)
}
\seealso{
\link{TCC-class}
}
\keyword{packages}
