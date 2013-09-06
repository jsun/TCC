\name{getNormalizedData}
\alias{getNormalizedData}
\title{Obtain normalized count data}
\description{
This function generates normalized count data from both original count
data and calculated normalization factors. 
}
\usage{getNormalizedData(tcc)}
\arguments{
    \item{tcc}{\link{TCC-class} object.}
}
\details{
This function is generally used after the \code{\link{calcNormFactors}}
function that calculates normalization factors. 
The normalized data is calculated using both the original count data 
stored in the \code{count} field and the normalization factors
stored in the \code{norm.factors} field in the \link{TCC-class} object.
}
\value{
A numeric matrix containing normalized count data.
}
\examples{
# Note that the hypoData has non-DEGs at 201-1000th rows
nonDEG <- 201:1000
data(hypoData)
summary(hypoData[nonDEG, ])
group <- c(1, 1, 1, 2, 2, 2)

# Obtaining normalized count data after performing the 
# DEGES/edgeR normalization method, i.e., DEGES/edgeR-normalized data
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
normalized.count <- getNormalizedData(tcc)
summary(normalized.count[nonDEG, ])

# Obtaining normalized count data after performing the 
# TMM normalization method (Robinson and Oshlack, 2010), 
# i.e., TMM-normalized data
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", iteration = 0)
normalized.count <- getNormalizedData(tcc)
summary(normalized.count[nonDEG, ])
}
\keyword{methods}
