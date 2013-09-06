\name{calcAUCValue}
\alias{calcAUCValue}
\title{Calculate AUC value from a TCC-class object}
\description{
This function calculates AUC (Area under the ROC curve) value from 
a \link{TCC-class} object for simulation study.
}
\usage{calcAUCValue(tcc)}
\arguments{
    \item{tcc}{\link{TCC-class} object having values in both \code{stat$rank}
               and \code{simulation$trueDEG} fields.}
}
\details{
This function is generally used after the \code{\link{estimateDE}} function 
that estimates \eqn{p}-values (and the derivatives such as the \eqn{q}-values 
and the ranks) for individual genes based on the
statistical model for differential expression (DE) analysis. 
In case of the simulation analysis, we know which genes are 
DEGs or non-DEGs in advance and the information is stored in 
the \code{simulation$trueDEG} field of the \code{\link{TCC-class}} 
object \code{tcc}  (i.e., \code{tcc$simulation$trueDEG}).
The \code{\link{calcAUCValue}} function calculates the AUC value 
between the ranked gene list obtained by
the \code{\link{estimateDE}} function and the truth 
obtained by the \code{\link{simulateReadCounts}} function. 
A well-ranked gene list should have a high AUC value 
(i.e., high sensitivity and specificity).
}
\value{
numeric scalar.
}
\examples{
# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) with biological replicates.
# the first 200 genes are DEGs, where 180 are up in G1.
# The DE analysis is performed by an exact test in edgeR coupled
# with the DEGES/edgeR normalization factors.
tcc <- simulateReadCounts(Ngene = 1000, PDEG = 0.2,
                         DEG.assign = c(0.9, 0.1),
                         DEG.foldchange = c(4, 4), 
                         replicates = c(3, 3))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
calcAUCValue(tcc)


# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) without replicates.
# the levels of DE are 3-fold in G1 and 7-fold in G2
# The DE analysis is performed by an negative binomial test in
# DESeq coupled with the DEGES/DESeq normalization factors.
tcc <- simulateReadCounts(Ngene = 1000, PDEG = 0.2,
                         DEG.assign = c(0.9, 0.1),
                         DEG.foldchange = c(3, 7), 
                         replicates = c(1, 1))
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "deseq", FDR = 0.1)
calcAUCValue(tcc)
}
\keyword{methods}
