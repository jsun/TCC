\name{ROKU}
\alias{ROKU}
\title{detect tissue-specific (or tissue-selective) patterns from
microarray data with many kinds of samples}
\description{
ROKU is a method for detecting tissue-specific (or tissue-selective)
patterns from gene expression data for many tissues (or samples).
ROKU (i) ranks genes according to their overall tissue-specificity
using Shannon entropy after data processing and (ii) detects tissues
specific to each gene if any exist using an Akaike's information
criterion (AIC) procedure.
}
\usage{
ROKU(data, upper.limit = 0.25, sort = FALSE)
}
\arguments{
    \item{data}{numeric matrix or data frame containing microarray data
          (on log2 scale), where each row indicates the gene or probeset ID, 
          each column indicates the tissue, and each cell indicates a
          (log2-transformed) expression value of the gene in the tissue.
          Numeric vector can also be accepted for a single gene expression
          vector.}
    \item{upper.limit}{numeric value (between 0 and 1) specifying the maximum
          percentage of tissues (or samples) as outliers to each gene.}
    \item{sort}{logical. If \code{TRUE}, results are sorted in descending
          order of the entropy scores.}
}
\details{
As shown in Figure 1 in the original study of ROKU (Kadota et al., 2006),
Shannon entropy \eqn{H} of a gene expression vector
(\eqn{x_{1}, x_{2}, ..., x_{N}}) for \eqn{N} tissues can range
from zero to \eqn{log_{2}N}, with the value 0 for genes expressed in a
single tissue and \eqn{log_{2}N} for genes expressed uniformly in all the
tissues. Researchers therefore rely on the low entropy score for the
identification of tissue-specific patterns.
However, direct calculation of the entropy for raw gene expression vector
works well only for detecting tissue-specific patterns when over-expressed
in a small number of tissues but unexpressed or slightly expressed in others:
The \eqn{H} scores of tissue-specific patterns such as
\eqn{(8,8,2,8,8,8,8,8,8,8)}
for the 3rd tissue-specific down-regulation (see the Figure 1e) are close
to the maximum value (\eqn{log_{2}N=3.32} when \eqn{N=10}) and cannot identify
such patterns as tissue-specific. To detect various kinds of tissue-specific
patterns by low entropy score, ROKU processes the original gene expression
vector and makes a new vector (\eqn{x_{1'}, x_{2'}, ..., x_{N'}}).
The data processing is done by subtracting the one-step Tukey biweight and
by taking the absolute value. In case of the above example, ROKU calculates
the \eqn{H} score from the processed vector \eqn{(0,0,6,0,0,0,0,0,0,0)},
giving very low score (from \eqn{H = 3.26} before processing to \eqn{H' = 0}
after processing). A major characteristic of ROKU is, therefore,
to be able to rank various tissue-specific patterns by using the modified
entropy scores.

Note that the modified entropy does not explain to which tissue a gene is
specific, only measuring the degree of overall tissue specificity of the gene.
ROKU employs an AIC-based outlier detection method (Ueda, 1996).
Consider, for example, a hypothetical mixed-type of tissue-selective expression
pattern \eqn{(1.2, 5.1, 5.2, 5.4, 5.7, 5.9, 6.0, 6.3, 8.5, 8.8)} where we
imagine a total of three tissues are specific (down-regulated in tissue1;
up-regulated in tissues 9 and 10). The method first normalize the expression
values by subtracting the mean and dividing by the standard deviation
(i.e., \eqn{z}-score transformation), then sorted in order of increasing
magnitude by\cr
\eqn{(-2.221, -0.342, -0.294, -0.198, -0.053, 0.043, 0.092, 0.236, 1.296,
1.441)}. The method evaluates various combinations of outlier candidates
starting from both sides of the values: model1 for non-outlier,
model2 for one outlier for high-side, model3 for two outliers for high-side,
..., model\eqn{x} for one outlier for down-side, ..., modely for two outliers for
both up- and down sides, and so on. Then, it calculates AIC-like statistic
(called \eqn{U}) for each combination of model and search the best combination
that achieves the lowest \eqn{U} value and is termed the minimum AIC estimate
(MAICE). Since the upper.limit value corresponds to the maximum number of the
outlier candidates, it decides the number of combinations. The AIC-based
method output a vector (1 for up-regulated outliers, -1 for down-regulated
outliers, and 0 for non-outliers) that corresponds to the input vector.
For example, the method outputs a vector \eqn{(-1, 0, 0, 0, 0, 0, 0, 0, 1, 1)}
when using \code{upper.limit = 0.5} and \eqn{(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0)}
when using \code{upper.limit = 0.25} (as default).
See the Kadota et al., 2007 for detailed discussion about the effect of
different parameter settings.
}
\value{
A list containing following fields:
    \item{outlier}{A numeric matrix when the input \code{data} are data frame
          or matrix. A numeric vector when the input \code{data} are numeric
          vector.
          Both matrix or vector consist of 1, -1, and 0: 1 for over-expressed
          outliers, -1 for under-expressed outliers, and 0 for non-outliers.}
    \item{modH}{A numeric vector when the input \code{data} are data frame or
          matrix. A numeric scalar when the input \code{data} are numeric
          vector. Both vector or scalar consist of modified entropy (\eqn{H'})
          score(s) calculated from a processed gene expression vector.}
    \item{H}{A numeric vector when the input \code{data} are data frame or
          matrix. A numeric scalar when the input \code{data} are numeric
          vector. Both vector or scalar consist of original entropy (\eqn{H})
          score(s) calculated from an original gene expression vector.}
    \item{rank}{A numeric vector or scalar consisting of the rank(s) of
          \code{modH}.}
    \item{Tbw}{a numeric vector or scalar consisting of one-step Tukey's
          biweight as an iteratively reweighted measure of central tendency.
          This value is in general similar to median value and the same as the
          output of \code{tukey.biweight} with default parameter settings
          in \code{affy} package. The data processing is done by
          subtracting this value for each gene expression vector and by
          taking the absolute value.}
}
\references{
Kadota K, Konishi T, Shimizu K: Evaluation of two outlier-detection-based
methods for detecting tissue-selective genes from microarray data.
Gene Regulation and Systems Biology 2007, 1: 9-15.

Kadota K, Ye J, Nakai Y, Terada T, Shimizu K:
ROKU: a novel method for identification of tissue-specific genes.
BMC Bioinformatics 2006, 7: 294.

Kadota K, Nishimura SI, Bono H, Nakamura S, Hayashizaki Y, Okazaki Y,
Takahashi K: Detection of genes with tissue-specific expression patterns
using Akaike's Information Criterion (AIC) procedure.
Physiol Genomics 2003, 12: 251-259.

Ueda T. Simple method for the detection of outliers.
Japanese J Appl Stat 1996, 25: 17-26.
}
\examples{
data(hypoData_ts)

result <- ROKU(hypoData_ts)
}
