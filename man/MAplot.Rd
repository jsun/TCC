\name{MAplot}
\alias{MAplot}
\title{plot a MA plot}
\description{
This method plots a MA plot based on the exactTestafterTbT.
This function will be obsoleted.
}
\usage{
MAplot(datalist, FDR_threshold = 0.01)
}
\arguments{
    \item{datalist}{The output from exacTestafterTbT}
    \item{FDR_threshold}{Points below the threshold will be plotted in red.}
}
\examples{
\dontrun{
sample <- NBsample()
out <- exactTestafterTbT(paste("gene",1:nrow(sample$data), sep=""), 
                        sample$data, c(1,1,1,2,2,2)) 
MAplot(out)
}
}
\seealso{\code{\link{edgeR}}}
