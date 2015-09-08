\name{makeFCMatrix}
\alias{makeFCMatrix}
\title{Generate the fold change matrix for simulating count data}
\description{
Generate the fold change matrix for simulating count data under the
Gammma distribution
}
\usage{
makeFCMatrix(Ngene = 10000, PDEG = 0.20, DEG.assign = NULL,
             replicates = NULL, fc.params = NULL)
}
\arguments{
    \item{Ngene}{numeric scalar specifying the number of genes.}
    \item{PDEG}{numeric scalar specifying the proportion of
          differentially expressed genes (DEGs).}
    \item{DEG.assign}{numeric vector specifying the proportion of DEGs up- or
          down-regulated in individual groups to be compared. The number of
          elements should be the same as that of \code{replicates} if
          \code{replicates} is specified. The indication of \code{replicates}
          means a single-factor experimental design.}
    \item{replicates}{numeric vector indicating the numbers of (biological)
          replicates for individual groups compared. Ignored if \code{group}
          is specified.}
    \item{fc.params}{foldchanges of DEGs are randomly sampled from
          \eqn{f + \Gamma(a, b)} where \eqn{a} is a shape and \eqn{b} is a scale
          of Gamma distribution. The default values are \eqn{f = 1.2},
          \eqn{a = 2}, and \eqn{b = 0.5}}
}
\details{
\code{makeFCMatrix} function is a function for generating the foldchanges of
DEGs. The foldchanges are randomly sampled from \eqn{f + \Gamma(a, b)} where
\eqn{a} is a shape and \eqn{b} is a scale of Gamma distribution.
}
\value{matrix}
\examples{
fc.matrix <- makeFCMatrix()
}
\keyword{methods}
