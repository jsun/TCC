\name{WAD}
\alias{WAD}
\title{Calculate WAD statistic for individual genes}
\description{
This function performs WAD method to identify differentially expressed genes
(DEGs) from two-group gene expression data. A high absolute value for the WAD
statistic is evident of a high degree of differential expression.
}
\usage{
WAD(data, group, log.scale = FALSE, floor.value = 1, sort = FALSE)
}
\arguments{
    \item{data}{numeric matrix or data frame containing count data or
          microarray data, where each row indicates the gene (or transcript
          or probeset ID), each column indicates the sample (or library),
          and each cell indicates the expression value (i.e., number of counts
          or signal intensity) of the gene in the sample.}
    \item{group}{numeric vector indicating the experimental group for each
          sample (or library).}
    \item{log.scale}{logical. If \code{TRUE}, the data are log scaled (base = 2)
          after performing the \code{floor.value} setting. The default is
          \code{log.scale = FALSE}.}
    \item{floor.value}{numeric scalar specifying a floor value for the data.
          The default is \code{floor.value = 1}, i.e., values less than 1 are
          replaced by 1. Ignored if \code{log.scale = FALSE}.}
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

wad <- WAD(nakai, group,
           log.scale = TRUE, floor.value = 1, sort = TRUE)
}
