\name{calcAUCValue}
\alias{calcAUCValue}
\title{Calculate AUC value from a TCC-class object}
\description{
This function calculates AUC (Area under the ROC curve) value from 
a \link{TCC-class} object for simulation study.
}
\usage{calcAUCValue(tcc, t = 1)}
\arguments{
    \item{tcc}{\link{TCC-class} object having values in both \code{stat$rank}
               and \code{simulation$trueDEG} fields.}
    \item{t}{numeric value (between 0 and 1) specifying the FPR (i.e., the
             \eqn{x}-axis of ROC curve). AUC value is calculated from 0 to
             \code{t}. The default is 1.}
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
# the first 200 genes are DEGs, where 180 are up-regulated in G1.
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

}
\keyword{methods}
