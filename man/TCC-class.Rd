\name{TCC-class}
\docType{class}
\alias{TCC-class}
\alias{show}
\alias{show,TCC-method}
\alias{names}
\alias{names,TCC-method}
\alias{length}
\alias{length,TCC-method}
\alias{[}
\alias{[,TCC-method}
\alias{subset,TCC-subset}
\title{A container for storing information used in TCC}
\description{
This is a container class for TCC. This class initially contains 
count data matrix and some information for the analysis of count data. 
It also provides further fields that are populated during the analysis.
}
\details{
This class is implemented as an R5 reference class.
Functions calling such methods copies the object prior to
calling the method to keep the semantics of functional programming.
This class can be created by the generic \code{new} function with 
count data and associated information of experimental design.

The values (defaults to all 1) in the \code{norm.factors} field 
will be changed after performing the \code{\link{calcNormFactors}} function.
The \code{DEGES} field stores information related to our DEGES-based
normalization pipeline after performing the \code{\link{calcNormFactors}}
function.
The \code{stat} and \code{estimatedDEG} fields store results after performing
the \code{\link{estimateDE}} function.
The \code{simulation} field stores parameters
used when performing the \code{\link{simulateReadCounts}} function.
}
\section{Fields}{
This class contains the following fields:
    \describe{
        \item{count}{numeric matrix containing count data. }
        \item{gene_id}{character vector indicating the identifier of the count
              unit, usually gene.}
        \item{group}{data frame indicating information about experimental
              design.}
        \item{norm.factors}{numeric vector containing normalization factors
              (default to 1).}
        \item{stat}{list for storing results after the execution of the 
              \code{\link{calcNormFactors}} and
              \code{\link{estimateDE}} functions.}
        \item{estimatedDEG}{numeric vector as a placeholder for indicating
              either DEGs (flagged as "1") or non-DEGs (as "0") for individual
              genes. The values in this field will be populated after
              the execution of the \code{\link{estimateDE}} function.}
        \item{simulation}{list. This field is only used for analyzing 
              simulation data.}
        \item{DEGES}{list for storing the information about normalization 
              steps.}
    }
}
\examples{
tcc <- simulateReadCounts(Ngene = 10000, PDEG = 0.2, DEG.assign = c(0.8, 0.2),
                          DEG.foldchange = c(4, 4), replicates = c(3, 3))

# Check the TCC class object.
tcc

# Check the fields of TCC class object.
names(tcc)
head(tcc$count)

# Check the normalization factors.
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors

# Check the p-values and q-values.
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
tcc

# Compare the breakdowns of estimated DEGs with the truth.
head(tcc$estimatedDEG)
head(tcc$simulation$trueDEG)

# M-A plotting.
plot(tcc)
}
\keyword{classes}
