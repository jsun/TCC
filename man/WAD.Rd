\name{WAD}
\alias{WAD}
\title{Calculate WAD statistic for individual genes}
\description{
This function performs WAD method to identify differentially expressed genes
(DEGs) from two-group gene expression data. A high absolute value for the WAD
statistic is evident of a high degree of differential expression.
}
\usage{
WAD(data, group, logged = FALSE, floor = 1, sort = FALSE)
}
\arguments{
    \item{data}{numeric matrix or data frame containing count data or
          microarray data, where each row indicates the gene (or transcript
          or probeset ID), each column indicates the sample (or library),
          and each cell indicates the expression value (i.e., number of counts
          or signal intensity) of the gene in the sample.}
    \item{group}{numeric vector indicating the experimental group for each
          sample (or library).}
    \item{logged}{logical. If \code{TRUE}, the input data are regarded as 
          log2-transformed. If \code{FALSE}, the log2-transformation is
          performed after the floor setting. The default is 
          \code{logged = FALSE}.}
    \item{floor}{numeric scalar (> 0) specifying the floor value for
          taking logarithm. The default is \code{floor = 1}, indicating that
          values less than 1 are replaced by 1. Ignored if
          \code{logged = TRUE}.}
    \item{sort}{logical. If \code{TRUE}, the retrieved results are sorted
          in order of the rank of absolute WAD statistic.
          If \code{FALSE}, the results are retrieved by the original order.}
}
\value{
A numeric vector of WAD statistic for individual genes
}
\references{
Kadota K, Nakai Y, Shimizu K: A weighted average difference method for
detecting differentially expressed genes from microarray data.
Algorithms Mol Biol. 2008, 3: 8.
}
\examples{
data(nakai)
group <- c(1, 1, 1, 1, 2, 2, 2, 2)

wad <- WAD(nakai, group, logged = TRUE, sort = TRUE)
}
