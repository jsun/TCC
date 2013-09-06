\name{hypoData_ts}
\docType{data}
\alias{hypoData_ts}
\title{A sample microarray data for detecting tissue-specific patterns.}
\description{
A hypothetical micoarray data consisting of eight rows (genes) and ten
columns (tissues). The expression patterns are quite similar to those
in figure 1 in Kadota et al., 2006.
}
\usage{data(hypoData_ts)}
\format{\code{hypoData_ts} is a matrix of dimension eight times ten.}
\details{
The \code{hypoData_ts} is designed for explaining the performance of
\code{\link{ROKU}} that identify tissue-specific expression patterns.
The \code{hypoData_ts} contains a total of eight genes having various
expression patterns across ten tissues: (1) 'up-type' genes selectively
over-expressed in a small number of tissues but unexpressed ("gene1"),
slightly expressed ("gene3"), and moderately expressed ("gene4"),
(2) 'down-type' genes selectively under-expressed ("gene5"),
and (3) 'mixed-type' genes selectively over- and under-expressed in
some tissues ("gene6"). The other genes are not tissue-specific genes
(i.e., "gene2", "gene7", and "gene8").
}
\references{
Kadota K, Ye J, Nakai Y, Terada T, Shimizu K:
ROKU: a novel method for identification of tissue-specific genes.
BMC Bioinformatics 2006, 7: 294.
}
\examples{
data(hypoData_ts)
}
\keyword{datasets}
