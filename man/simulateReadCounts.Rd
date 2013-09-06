\name{simulateReadCounts}
\alias{simulateReadCounts}
\title{Generate simulation data from negative binomial (NB) distribution}
\description{
This function generates simulation data with arbitrary defined
experimental condition.
}
\usage{
simulateReadCounts(Ngene = 10000, PDEG = 0.20, DEG.assign = NULL,
                   DEG.foldchange = NULL, replicates = NULL, group = NULL)
}
\arguments{
    \item{Ngene}{numeric scalar specifying the number of genes.}
    \item{PDEG}{numeric scalar specifying the proportion of
          differentially expressed genes (DEGs).}
    \item{DEG.assign}{numeric vector specifying the proportion of DEGs up- or
          down-regulated in individual groups to be compared. The number of
          elements should be the same as that of \code{replicates} if
          \code{replicates} is specified. The indication of \code{replicates}
          means a single-factor experimental design. The number of elements in
          \code{DEG.assign} should be the same as the number of columns in
          \code{DEG.foldchange}. Both \code{DEG.foldchange} as data frame
          and \code{group} should simultaneously be specified and those
          indication means a multi-factor experimental design.}
    \item{DEG.foldchange}{numeric vector for single-factor experimental design
          and data frame for multi-factor experimental design. Both
          \code{DEG.foldchange} as numeric vector and \code{replicates} should
          simultaneously be specified for single-factor experimental design.
          The \eqn{i}-th element in \code{DEG.foldchange} vector indicates the
          degree of fold-change for Group \eqn{i}. The default is
          \code{DEG.foldchange = c(4, 4)}, indicating that the levels of DE
          are four-fold in both groups.\cr
          Both \code{DEG.foldchange} as data frame and \code{group} should
          simultaneously be specified for multi-factor experimental design.
          Numeric values in the \code{DEG.foldchange} object indicate the
          degree of fold-change for individual conditions or factors.}
    \item{replicates}{numeric vector indicating the numbers of (biological)
          replicates for individual groups compared. Ignored if \code{group}
          is specified.}
    \item{group}{data frame specifying the multi-factor experimental design.}
}
\details{
The empirical distribution of read counts 
used in this function is calculated from a RNA-seq dataset 
obtained from \emph{Arabidopsis} data 
(three biological replicates for both the treated and non-treated samples), 
the \code{arab} object, in NBPSeq package (Di et al., 2011). 
The overall design about the simulation conditions introduced 
can be viewed as a pseudo-color image by the 
\code{\link{plotFCPseudocolor}} function.
}
\value{
A \link{TCC-class} object containing following fields:
    \item{count}{numeric matrix of simulated count data.}
    \item{group}{data frame indicating which group (or condition or factor)
          each sample belongs to.}
    \item{norm.factors}{numeric vector as a placeholder for 
                        normalization factors.}
    \item{stat}{list for storing results after the execution of
                the \code{\link{calcNormFactors}} 
                (and \code{\link{estimateDE}}) function.}
    \item{estimatedDEG}{numeric vector as a placeholder for indicating
                        which genes are up-regulated in particular group 
                        compared to the others. The values in this field
                        will be populated after the execution of the 
                        \code{\link{estimateDE}} function.}
    \item{simulation}{list containing four fields: \code{trueDEG},
                      \code{DEG.foldchange}, \code{PDEG}, and \code{params}. 
                      The \code{trueDEG} field (numeric vector) stores 
                      information about DEGs: 0 for non-DEG, 1 for 
                      DEG up-regulated in Group 1, 2 for DEG up-regulated 
                      in Group 2, and so on. The information for
                      the remaining three fields is the same as those 
                      indicated in the corresponding arguments.}
}
\examples{
# Generating a simulation data for comparing two groups
# (G1 vs. G2) without replicates (single-factor experimental design). 
# the levels of DE are 3-fold in G1 and 7-fold in G2
tcc <- simulateReadCounts(Ngene = 10000, PDEG = 0.2, 
                         DEG.assign = c(0.9, 0.1),
                         DEG.foldchange = c(3, 7),
                         replicates = c(1, 1))
dim(tcc$count)
head(tcc$count)
str(tcc$simulation)
head(tcc$simulation$trueDEG)


# Generating a simulation data for comparing three groups
# (G1 vs. G2 vs. G3) with biological replicates
# (single-factor experimental design).
# the first 3000 genes are DEGs, where the 70%, 20%, and 10% are
# up-regulated in G1, G2, G3, respectively. The levels of DE are
# 3-, 10-, and 6-fold in individual groups.
tcc <- simulateReadCounts(Ngene = 10000, PDEG = 0.3, 
                         DEG.assign = c(0.7, 0.2, 0.1),
                         DEG.foldchange = c(3, 10, 6), 
                         replicates = c(2, 4, 3))
dim(tcc$count)
head(tcc$count)
str(tcc$simulation)
head(tcc$simulation$trueDEG)


# Generating a simulation data consisting of 10,000 rows (i.e., Ngene = 10000)
# and 8 columns (samples) for two-factor experimental design
# (condition and time). The first 3,000 genes are DEGs (i.e., PDEG = 0.3).
# Of the 3,000 DEGs, 40% are differentially expressed in condition (or GROUP) "A"
# compared to the other condition (i.e., condition "B"), 40% are differentially
# expressed in condition (or GROUP) "B" compared to the other condition
# (i.e., condition "A"), and the remaining 20% are differentially expressed at
# "10h" in association with the second factor: DEG.assign = c(0.4, 0.4, 0.2).
# The levels of fold-change are (i) 2-fold up-regulation in condition "A" for
# the first 40% of DEGs, (ii) 4-fold up-regulation in condition "B" for the
# second 40%, and (iii) 0.4- and 0.6-fold up-regulation at "10h" in "A" and
# 5-fold up-regulation at "10h" in "B".

group <- data.frame(
   GROUP = c("A", "A", "A", "A", "B", "B", "B", "B"),
   TIME  = c("2h", "2h", "10h", "10h", "2h", "2h", "10h", "10h")
)
DEG.foldchange <- data.frame(
   FACTOR1.1 = c(2, 2, 2, 2, 1, 1, 1, 1),
   FACTOR1.2 = c(1, 1, 1, 1, 4, 4, 4, 4),
   FACTOR2   = c(1, 1, 0.4, 0.6, 1, 1, 5, 5)
)
tcc <- simulateReadCounts(Ngene = 10000, PDEG = 0.3,
                          DEG.assign = c(0.4, 0.4, 0.2),
                          DEG.foldchange = DEG.foldchange,
                          group = group)
tcc
}
\keyword{methods}
