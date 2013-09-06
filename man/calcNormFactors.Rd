\name{calcNormFactors}
\docType{methods}
\alias{calcNormFactors,TCC-method}
\alias{calcNormFactors,DGEList-method}
\alias{calcNormFactors}
\title{Calculate normalization factors}
\description{
This function calculates normalization factors using a specified 
multi-step normalization method from a \link{TCC-class} object.
The procedure can generally be described as the 
\eqn{STEP1-(STEP2-STEP3)_{n}}{STEP1-(STEP2-STEP3)n} pipeline.
}
\usage{
\S4method{calcNormFactors}{TCC}(tcc, norm.method = NULL, test.method = NULL,
                iteration = TRUE,  FDR = NULL, floorPDEG = 0.05, 
                increment = FALSE, ...)
}
\arguments{
    \item{tcc}{\link{TCC-class} object.}
    \item{norm.method}{character specifying a normalization method used in
          both the STEP1 and STEP3. Possible values are \code{"tmm"} for the
          TMM normalization method implemented in the edgeR package,
          \code{"edger"} (same as \code{"tmm"}), and \code{"deseq"} for the
          method implemented in the DESeq package. The default is \code{"tmm"}
          when analyzing the count data with multiple replicates
          (i.e., \code{min(table(tcc$group[, 1])) > 1}) 
          and \code{"deseq"} when analyzing the count data without replicates
          \cr(i.e., \code{min(table(tcc$group[, 1])) == 1}).}
    \item{test.method}{character specifying a method for identifying
          differentially expressed genes (DEGs) used in STEP2:
          one of \code{"edger"}, \code{"deseq"}, \code{"bayseq"},
          \code{"ebseq"}, \code{"samseq"}, and \code{"wad"}.
          See the "Details" filed in \code{\link{estimateDE}} for detail.
          The default is \code{"edger"} when analyzing the count data with
          multiple replicates (i.e., \code{min(table(tcc$group[, 1])) > 1}) 
          and \code{"deseq"} when analyzing the count data without replicates 
          (i.e., \code{min(table(tcc$group[, 1])) == 1))}.}
    \item{iteration}{logical or numeric value specifying the number of
          iteration (\eqn{n}) in the proposed normalization pipeline: the
          \eqn{STEP1-(STEP2-STEP3)_{n}}{STEP1-(STEP2-STEP3)n} pipeline.
          If \code{FALSE} or \code{0} is specified, the normalization pipeline
          is performed only by the method in STEP1. If \code{TRUE} or \code{1}
          is specified, the three-step normalization pipeline is performed.
          Integers higher than \code{1} indicate the number of iteration in
          the pipeline.}
    \item{FDR}{numeric value (between 0 and 1) specifying the threshold for
          determining potential DEGs after STEP2.}
    \item{floorPDEG}{numeric value (between 0 and 1) specifying the minimum
          value to be eliminated as potential DEGs before performing STEP3.}
    \item{increment}{logical value. if \code{increment = TRUE}, the DEGES
          pipeline will perform again from the current iterated result.}
    \item{...}{arguments to identify potential DEGs at STEP2. See the
          "Arguments" field in \code{\link{estimateDE}} for details.}
}
\details{
The \code{\link{calcNormFactors}} function is the main function in the
TCC package. 
Since this pipeline employs the DEG identification method at STEP2, 
our multi-step strategy can eliminate the negative effect of potential DEGs
before the second normalization at STEP3.
To fully utilize the DEG elimination strategy (DEGES), we strongly recommend
not to use \code{iteration = 0} or \code{iteration = FALSE}. 
This function internally calls functions implemented in the other R packages
according to the specified value. 

\itemize{
    \item \code{norm.method = "tmm"}\cr
          The \code{\link[edgeR]{calcNormFactors}} function implemented
          in edgeR is used for obtaining the TMM normalization factors
          at both STEP1 and STEP3.
    \item \code{norm.method = "deseq"}\cr
          The \code{\link[DESeq]{estimateSizeFactors}} function 
          implemented in DESeq is used for obetaining the size factors
          at both STEP1 and STEP3.
          The size factors are internally converted to normalization factors
          that are comparable to the TMM normalization factors.
}
}
\value{
After performing the \code{calcNormFactors} function,
the calculated normalization factors are populated in the
\code{norm.factors} field (i.e., \code{tcc$norm.factors}).
Parameters used for DEGES normalization (e.g., potential DEGs
identified in STEP2, execution times for the identification, etc.)
are stored in the DEGES field (i.e., \code{tcc$DEGES}) as follows:

    \item{iteration}{the iteration number \eqn{n} for
          \bold{the STEP1 - (STEP2 - STEP3)\eqn{_{n}} pipeline}.}
    \item{pipeline}{the DEGES normalization pipeline.}
    \item{threshold}{it stores
          (i) the type of threshold (\code{threshold$type}),
          (ii) the threshold value (\code{threshold$input}),
          and (iii) the percentage of potential DEGs actually
          used (\code{threshold$PDEG}).
          These values depend on whether the percentage
          of DEGs identified in STEP2 is higher or lower to the value
          indicated by \code{floorPDEG}.
          Consider, for example, the execution of \code{calcNormFactors}
          function with "\code{FDR = 0.1} and \code{floorPDEG = 0.05}".
          If the percentage of DEGs identified in STEP2 satisfying 
          \code{FDR = 0.1} was 0.14
          (i.e., higher than the \code{floorPDEG} of 0.05),
          the values in the threshold fields will be
          \code{threshold$type = "FDR"}, \code{threshold$input = 0.1},
          and \code{threshold$PDEG = 0.14}.
          If the percentage (= 0.03) was lower than the predefined
          \code{floorPDEG} value of 0.05, the values in the threshold fields
          will be \code{threshold$type = "floorPDEG"}, 
          \code{threshold$input = 0.05}, and \code{threshold$PDEG = 0.05}.}
    \item{potDEG}{numeric binary vector (0 for non-DEG or 1 for DEG)
          after the evaluation of the percentage of DEGs identified in STEP2 with
          the predefined \code{floorPDEG} value. If the percentage
          (e.g., 2\%) is lower than the \code{floorPDEG} value (e.g., 17\%),
          17\% of elements become 1 as DEG.}
    \item{prePotDEG}{numeric binary vector
          (0 for non-DEG or 1 for DEG) before the evaluation of the percentage
          of DEGs identified in STEP2 with the predefined \code{floorPDEG}
          value. Regardless of the \code{floorPDEG} value, the percentage of
          elements with 1 is always the same as that of DEGs identified in
          STEP2.}
    \item{execution.time}{computation time required for normalization.}
}
\examples{
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)

# calculating normalization factors using the DEGES/edgeR method 
# (the TMM-edgeR-TMM pipeline)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors

# calculating normalization factors using the iterative DEGES/edgeR method 
# (iDEGES/edgeR) with n = 3
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 3, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors

# calculating normalization factors for simulation data without replicates
tcc <- simulateReadCounts(replicates = c(1, 1))
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors
}
