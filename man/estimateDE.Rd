\name{estimateDE}
\alias{estimateDE}
\title{Estimate degrees of differential expression (DE) for individual genes}
\usage{
estimateDE(tcc, test.method, FDR, 
           dispersion, fit0, fit1, design,
           contrast, coef, comparison, samplesize, floor.value = 1, cl)
}
\description{
This function calculates \eqn{p}-values (or the related statistics) for
identifying differentially expressed genes (DEGs) from a
\link{TCC-class} object. 
\code{estimateDE} internally calls a specified method 
implemented in other R packages.
}
\arguments{
    \item{tcc}{\link{TCC-class} object.}
    \item{test.method}{character string specifying a method for identifying
          DEGs: one of \code{"edger"}, \code{"deseq"}, \code{"bayseq"},
          \code{"ebseq"}, \code{"samseq"}, and \code{"wad"}.
          See the "Details" field for detail.
          The default is \code{"edger"} when analyzing the count data with
          replicates (i.e., \code{min(table(tcc$group[, 1])) > 1}), 
          \code{"deseq"} when analyzing the count data without replicates 
          (i.e., \code{min(table(tcc$group[, 1])) == 1)}).}
          %, and \code{"bayseq"} when the \code{paired = TRUE} is specified.}
    \item{FDR}{numeric value (between 0 and 1) specifying the threshold
          for determining DEGs.}
    %\item{paired}{logical. If \code{TRUE}, the data is treated as paired data.}
    \item{dispersion}{numeric vector giving the dispersion of all genes 
          for analyzing the count data without replicates when  the
          \code{test.method = "edger"} is specified.
          See the \code{\link[edgeR]{exactTest}} function in edgeR for details.}
    \item{design}{numeric matrix giving the design matrix for the linear model.
          See the \code{\link[edgeR]{glmFit}} function in edgeR for details.}
    \item{contrast}{numeric vector specifying a contrast of the linear model
          coefficients to be tested equal to zero.
          See the \code{\link[edgeR]{glmLRT}} function in edgeR for details.}
    \item{coef}{integer or character vector indicating which coefficients
          of the linear model are to be tested equal to zero.
          See the \code{\link[edgeR]{glmLRT}} function in edgeR for details.}
    \item{fit0}{a formula for creating reduced model described in DESeq. 
          The left hand side must be \code{count}, the right hand side can
          involve any column of \code{tcc$group} is used as the model frame.
          See the \code{\link[DESeq]{fitNbinomGLMs}} function for details.}
    \item{fit1}{a formula for creating full model described in DESeq.
          The left hand side must be \code{count}, the right hand side can
          involve any column of \code{tcc$group} is used as the model frame.
          See the \code{\link[DESeq]{fitNbinomGLMs}} function for details.}
    \item{comparison}{numeric or character string identifying the columns in
          the \code{tcc$group[, 1]} for analysis. See the \code{group} argument
          in \code{\link[baySeq]{topCounts}} for details.}
    \item{samplesize}{integer specifying (i) the sample size for estimating the
          prior parameters if \code{test.method = "bayseq"} (defaults to 10000),
          (ii) the iteration number used in EM algorithm in EBSeq if 
          \code{test.method = "ebseq"} (defaults to 5),
          and (iii) the number of permutation in samr if
          \code{test.method = "samseq"} (defaults to 10).}
    \item{floor.value}{numeric scalar (>= 0) to replace the values higher
          than the \code{floor.value} if \code{test.method = "wad"}
          (defaults to 1).}
    \item{cl}{\code{snow} object when using multi processors if
          \code{test.method = "bayseq"} is specified.
          See the \code{\link[baySeq]{getPriors.NB}} function in baySeq
          for details.}
}
\details{
\code{estimaetDE} function is generally used after performing the
\code{\link{calcNormFactors}} function that calculates normalization factors.
\code{estimateDE} constructs a statistical model for differential expression
(DE) analysis with the calculated normalization factors and returns the
\eqn{p}-values (or the derivatives). The individual functions in the other
packages are internally called according to the specified
\code{test.method} parameter.

\itemize{
    \item \code{test.method = "edger"}\cr
          There are two approaches (i.e., exact test and GLM) to identify DEGs
          in edgeR. The two approches are implmented in TCC. As a default,
          the exact test approach is used for two-group data, 
          and GLM approach is used for multi-group or multi-factor data. 
          However, if \code{design} and the one of \code{coef} or 
          \code{contrast} are given, the GLM approach will be used for
          two-group data. Currently, edgeR package does not support data
          without replicates.\cr
          If the exact test approach is used, 
          \code{\link[edgeR]{estimateCommonDisp}},
          \code{\link[edgeR]{estimateTagwiseDisp}}, and
          \code{\link[edgeR]{exactTest}} are internally called.\cr
          If the GLM approach is used,
          \code{\link[edgeR]{estimateGLMCommonDisp}},
          \code{\link[edgeR]{estimateGLMTrendedDisp}},
          \code{\link[edgeR]{estimateGLMTagwiseDisp}},
          \code{\link[edgeR]{glmFit}}, and
          \code{\link[edgeR]{glmLRT}}
          are internally called.
    \item \code{test.method = "deseq"}\cr
          DESeq supports two approach (i.e. an exact test and 
          GLM approach) for identifying DEGs. As a default,
          the exact test is used for two-group data,
          and GLM approach is used for multi-group or multi-factor data.
          However, if \code{fit0} and \code{fit1} are given, the GLM approach
          will be used for two-group data.\cr
          If the exact test is used, 
          \code{\link[DESeq]{estimateDispersions}} and
          \code{\link[DESeq]{nbinomTest}} are internally called.\cr
          If the GLM approach is used, 
          \code{\link[DESeq]{estimateDispersions}},
          \code{\link[DESeq]{fitNbinomGLMs}}, and
          \code{\link[DESeq]{nbinomGLMTest}}
          are internally called.
    \item \code{test.method = "bayseq"}\cr
          \code{\link[baySeq]{getPriors.NB}} and
          \code{\link[baySeq]{getLikelihoods.NB}} in baySeq are internally
          called for identifying DEGs.
          %when \code{paired = FALSE}.
          %\code{\link[baySeq]{getPriors.BB}} and
          %\code{\link[baySeq]{getLikelihoods.BB}} in baySeq are called
          %when \code{paired = TRUE}.
          Since baySeq outputs posterior likelihoods instead of
          \eqn{p}-values, the \eqn{(1 - likelihood)} values are
          stored in the corersponding field
          (i.e., \code{tcc$stat$p.value}).
    \item \code{test.method = "ebseq"}\cr
          \code{\link[EBSeq]{EBTest}} and \code{\link[EBSeq]{EBMultiTest}}
          are internally called for two-group and multi-group data with
          replicates, respectively. Currently, both functions does not
          support data without replicates.
          Since EBSeq only outputs \eqn{q}-value,
          the all values in \code{tcc$stat$p.value} are \code{NA}.
    \item \code{test.method = "samseq"}\cr
          \code{\link[samr]{SAMseq}} with
          \code{resp.type = "Two class unpaired"} arugment
          in samr package is called to identify DEGs for two-group data,
          and \code{resp.type = "Multiclass"} for multi-group data.
          Since \code{\link[samr]{SAMseq}} outputs test statistics
          instead of \eqn{p}-values, 
          the \code{tcc$stat$p.value} and \code{tcc$stat$q.value}
          are  \code{NA}.
          Alternatively, the test statistics are stored in 
          \code{tcc$stat$testStat} field.
    \item \code{test.method = "wad"}\cr
          The \code{\link{WAD}} implemented in TCC is used for identifying
          DEGs. Since \code{\link{WAD}} outputs test statistics instead of
          \eqn{p}-values, the \code{tcc$stat$p.value} and
          \code{tcc$stat$q.value} are \code{NA}.
          Alternatively, the test statistics are stored in
          \code{tcc$stat$testStat} field.
}
}
\value{
A \code{\link{TCC-class}} object containing following fields:
    \item{stat$p.value}{numeric vector of \eqn{p}-values.}
    \item{stat$q.value}{numeric vector of \eqn{q}-values calculated
          based on the \eqn{p}-values using the \code{p.adjust} function
          with default parameter settings.}
    \item{stat$testStat}{numeric vector of test statistics if
          \code{test.method = "samseq"} or \code{"wad"} is specified.}
    \item{stat$rank}{gene rank in order of the \eqn{p}-values or
          test statistics.}
    \item{estimatedDEG}{numeric vector consisting of 0 or 1
                        depending on whether each gene is classified 
                        as non-DEG or DEG. The threshold for classifying
                        DEGs or non-DEGs is preliminarily given as the 
                        \code{FDR} argument.}
}
\examples{
# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) with biological replicates
# The DE analysis is performed by an exact test in edgeR coupled
# with the DEGES/edgeR normalization factors.
# For retrieving the summaries of DE results, we recommend to use
# the getResult function.
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
head(tcc$stat$p.value)
head(tcc$stat$q.value)
head(tcc$estimatedDEG)
result <- getResult(tcc)


# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) without replicates
# The DE analysis is performed by an negative binomial (NB) test
# in DESeq coupled with the DEGES/DESeq normalization factors.
data(hypoData)
group <- c(1, 2)
tcc <- new("TCC", hypoData[, c(1, 4)], group)
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "deseq", FDR = 0.1)
}
\keyword{methods}
