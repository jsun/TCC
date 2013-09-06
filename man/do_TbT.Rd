\name{do_TbT}
\alias{do_TbT}
\title{Calculate normalization factors for raw tag count data using a
multi-step normalization strategy called "TbT"}
\description{
This method performs TMM normalization, baySeq, and again TMM 
normalization to infer a good normalization factor as
described in Kadota et al. 2012, Algorithms Mol Biol 7:5.
This function will be obsoleted. 
Please use \link{TCC-class} based methods.
}
\usage{
do_TbT(data, data.cl, sample_num = 10000)
}
\arguments{
    \item{data}{The data matrix to be analysed. Numerical data only}
    \item{data.cl}{A vector describing the data class for columns in data}
    \item{sample_num}{Sample number for baysian estimation}
}
\examples{
\dontrun{
sample <- NBsample()
out <- do_TbT(sample$data, c(1, 1, 1, 2, 2, 2))
}
}
\seealso{\code{\link{edgeR}}}
