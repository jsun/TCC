\name{getResult}
\alias{getResult}
\title{Obtain the summaries of results after the differential
 expression analysis}
\description{
This function is generally used after the
\code{\link{estimateDE}} function. It retrieves the summaries of
differential expression (DE) results from \link{TCC-class} object.
The retrieved information includes \eqn{p}-values, \eqn{q}-values,
coordinates of M-A plot (i.e., M and A values), and so on.
}
\usage{getResult(tcc, sort = FALSE, floor = 0)}
\arguments{
    \item{tcc}{\link{TCC-class} object}
    \item{sort}{logical. If \code{TRUE}, the retrieved results are
                sorted in order of the \code{stat$rank} field in the
                \link{TCC-class} object. If \code{FALSE}, the results
                are retrieved by the original order.}
     \item{floor}{numeric scalar specifying a threshold for adjusting
                  low count data.}
}
\value{
A data frame object containing following fields:
    \item{gene_id}{character vector indicating the id of the count unit,
                   usually gene.}
    \item{a.value}{numeric vector of average expression level on log2
                   scale (i.e., A-value) for each gene across the compared
                   two groups. It corresponds to the \eqn{x} coordinate 
                   in the M-A plot.}
    \item{m.value}{numeric vector of fold-change on \eqn{\log_2} scale
                   (i.e., M-value) for each gene between the two groups
                   compared. It corresponds to the \eqn{y} coordinate in
                   the M-A plot.}
    \item{p.value}{numeric vector of \eqn{p}-values.}
    \item{q.value}{numeric vector of \eqn{q}-values calculated
                   based on the \eqn{p}-values using the \code{p.adjust}
                   function with default parameter settings.}
    \item{rank}{numeric vector of gene rank in order of the \eqn{p}-values.}
    \item{estimatedDEG}{numeric vector consisting of 0 or 1
                        depending on whether each gene is classified 
                        as non-DEG or DEG. The threshold for classifying 
                        DEGs or non-DEGs is preliminarily 
                        given when performing \code{\link{estimateDE}}.}
}
\examples{
# Obtaining DE results by an exact test in edgeR coupled with
# the DEGES/edgeR normalization factors
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)
head(result)

# Obtaining DE results by an negative binomial test in DESeq
# coupled with the iterative DEGES/DESeq normalization method
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "deseq", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)
head(result)
}
\keyword{methods}
