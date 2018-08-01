\name{estimateDE}
\alias{estimateDE}
\title{Estimate degrees of differential expression (DE) for individual genes}
\usage{
estimateDE(tcc, test.method, FDR, paired,
           full, reduced,                   # for DESeq, DESeq2
           design, contrast,                # for edgeR, DESeq2, voom
           coef,                            # for edgeR, voom
           group, cl,                       # for baySeq
           samplesize,                      # for baySeq, SAMseq
           logged, floor,                   # for WAD
           ...
)
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
          DEGs: one of \code{"edger"}, \code{"deseq"}, \code{"deseq2"},
          \code{"bayseq"}, \code{"voom"}, and \code{"wad"}.
          See the "Details" field for detail.
          The default is \code{"edger"} when analyzing the count data with
          replicates (i.e., \code{min(table(tcc$group[, 1])) > 1}), and
          \code{"deseq"} (2 group) and \code{"deseq2"} (more than 2 group)
          when analyzing the count data without replicates 
          (i.e., \code{min(table(tcc$group[, 1])) == 1}).}

    \item{FDR}{numeric value (between 0 and 1) specifying the threshold
          for determining DEGs.}

    \item{paired}{logical. If \code{TRUE}, the input data are regarded as
          (two-group) paired samples. If \code{FALSE}, the input data are
          regarded as unpaired samples. The default is \code{FALSE}.}

    \item{full}{a formula for creating full model described in DESeq and
          DESeq2. The right hand side can involve any column of \code{tcc$group}
          is used as the model frame.
          See the \code{\link[DESeq]{fitNbinomGLMs}} function in DESeq 
          for details, or \code{\link[DESeq2]{nbinomLRT}} function in DESeq2.}

    \item{reduced}{a formula for creating reduced model described in DESeq.
          The right hand side can involve any column of \code{tcc$group} is
          used as the model frame.
          See the \code{\link[DESeq]{fitNbinomGLMs}} function in DESeq 
          for details, or \code{\link[DESeq2]{nbinomLRT}} function in DESeq2.}

    \item{design}{the argument is used in edgeR, voom (limma) and DESeq2.
          For edgeR and voom, it should be the numeric matrix giving the
          design matrix for the generalized linear model.
          See the \code{\link[edgeR]{glmFit}} function in edgeR or
          the \code{\link[limma]{lmFit}} function in limma for details.
          For DESeq2, it should be a formula specifying the design of the
          experiment. See the \code{\link[DESeq2]{DESeqDataSet}} function
          in DESeq2 for details.}

    \item{contrast}{the argument is used in edgeR and DESeq2.
          For edgeR, numeric vector specifying a contrast of the linear model
          coefficients to be tested equal to zero.
          See the \code{\link[edgeR]{glmLRT}} function in edgeR for details.
          For DESeq2, the argument is same to \code{contrast} which used in
          DESeq2 package to retrive the results from Wald test. See the
          \code{\link[DESeq2]{results}} function in DESeq2 for details.}

    \item{coef}{integer or character vector indicating which coefficients
          of the linear model are to be tested equal to zero.
          See the \code{\link[edgeR]{glmLRT}} function in edgeR for details.}

    \item{group}{numeric or character string identifying the columns in
          the \code{tcc$group} for analysis. See the \code{group} argument
          of \code{\link[baySeq]{topCounts}} function in baySeq for details.}

    \item{cl}{\code{snow} object when using multi processors if
          \code{test.method = "bayseq"} is specified.
          See the \code{\link[baySeq]{getPriors.NB}} function in baySeq
          for details.}

    \item{samplesize}{integer specifying the sample size for estimating the
          prior parameters if \code{test.method = "bayseq"} (defaults to 10000).}

    \item{logged}{logical. If \code{TRUE}, the input data are regarded as 
          log2-transformed. If \code{FALSE}, the log2-transformation is
          performed after the floor setting. The default is 
          \code{logged = FALSE}.
          Ignored if \code{test.method} is not \code{"wad"}.}

    \item{floor}{numeric scalar (> 0) specifying the floor value for
          taking logarithm. The default is \code{floor = 1}, indicating that
          values less than 1 are replaced by 1. Ignored if
          \code{logged = TRUE}.
          Ignored if \code{test.method} is not \code{"wad"}.}

    \item{...}{further paramenters.}
}
\details{
\code{estimaetDE} function is generally used after performing the
\code{\link{calcNormFactors}} function that calculates normalization factors.
\code{estimateDE} constructs a statistical model for differential expression
(DE) analysis with the calculated normalization factors and returns the
\eqn{p}-values (or the derivatives). The individual functions in other
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
          two-group data. \cr
          If the exact test approach is used, 
          \code{\link[edgeR]{estimateCommonDisp}},
          \code{\link[edgeR]{estimateTagwiseDisp}}, and
          \code{\link[edgeR]{exactTest}} are internally called.\cr
          If the GLM approach is used,
          \code{\link[edgeR]{estimateGLMCommonDisp}},
          \code{\link[edgeR]{estimateGLMTrendedDisp}},\cr
          \code{\link[edgeR]{estimateGLMTagwiseDisp}},
          \code{\link[edgeR]{glmFit}}, and
          \code{\link[edgeR]{glmLRT}}
          are internally called.
    \item \code{test.method = "deseq"}\cr
          DESeq supports two approach (i.e. an exact test and 
          GLM approach) for identifying DEGs. As a default,
          the exact test is used for two-group data,
          and GLM approach is used for multi-group or multi-factor data.
          However, if \code{full} and \code{reduced} are given, the GLM approach
          will be used for two-group data.\cr
          If the exact test is used, 
          \code{\link[DESeq]{estimateDispersions}} and
          \code{\link[DESeq]{nbinomTest}} are internally called.\cr
          If the GLM approach is used, 
          \code{\link[DESeq]{estimateDispersions}},
          \code{\link[DESeq]{fitNbinomGLMs}}, and
          \code{\link[DESeq]{nbinomGLMTest}}
          are internally called.
    \item \code{test.method = "deseq2"}\cr
          \code{\link[DESeq2]{estimateDispersions}}, and 
          \code{\link[DESeq2]{nbinomWaldTest}} are internally called for 
          identifying DEGs.
          However, if \code{full} and \code{reduced} are given,
          the \code{\link[DESeq2]{nbinomLRT}} will be used.
    \item \code{test.method = "bayseq"}\cr
          \code{\link[baySeq]{getPriors.NB}} and
          \code{\link[baySeq]{getLikelihoods}} in baySeq are internally
          called for identifying DEGs.
          If \code{paired = TRUE}, 
          \code{\link[baySeq]{getPriors}} and
          \code{\link[baySeq]{getLikelihoods}} in baySeq are used.
    \item \code{test.method = "voom"}\cr
          \code{\link[limma]{voom}}, \code{\link[limma]{lmFit}}, and
          \code{\link[limma]{eBayes}} in limma are internally called
          for identifying DEGs. 
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
          \code{"wad"} is specified.}
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
