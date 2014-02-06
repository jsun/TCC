\name{clusterSample}
\alias{clusterSample}
\title{Perform hierarchical clustering for samples from expression data}
\description{
This function performs hierarchical clustering for samples
(tissues or columns) from expression data.
}
\usage{
clusterSample(data, dist.method = "spearman", hclust.method = "average",
              unique.pattern = TRUE)
}
\arguments{
    \item{data}{numeric matrix or data frame containing expression data
         (count data or microarray data), where each row indicates the
         gene (or transcript or probeset ID), each column indicates the
         sample (or library), and each cell indicates the expression value
         (i.e., number of counts or signal intensity) of the gene in the
         sample.}
    \item{dist.method}{character string specifying a type for correlation
         coefficient (\code{"spearman"} or \code{"pearson"}) used as
         distance. The default is \code{"spearman"}. The hierarchical
         clustering is performed using the distance
         (i.e., 1 - \code{"spearman"} correlation coefficient, by default).}
    \item{hclust.method}{character string specifying an agglomeration method
         used in \code{\link[stats]{hclust}} function: \code{"ward"},
         \code{"single"}, \code{"complete"}, \code{"average"},
         \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
         The default is \code{"average"}.}
    \item{unique.pattern}{logical. If \code{FALSE}, the input expression
         data are directly applied for clustering. If \code{TRUE} (default),
         the input data only having unique expression patterns are applied.)}
}
\value{
An object of class \code{\link[stats]{hclust}} which describes the tree
produced by the clustering process. See \code{\link[stats]{hclust}} for
details.
}
\examples{
# 1. Perform sample clustering with default options.
data(hypoData)
hc <- clusterSample(hypoData)
plot(hc)

# 2. Obtain the same result using the 'unique.pattern = FALSE' option.
data(hypoData)
keep <- as.logical(rowSums(hypoData) > 0)
data <- unique(hypoData[keep, ])
hc <- clusterSample(data, unique.pattern = FALSE)
plot(hc)
}
