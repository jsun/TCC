\name{arab}
\docType{data}
\alias{arab}
\title{Arabidopsis RNA-Seq data set}
\description{
This dataset was imported from NBPSeq package and the following
explanation is verbatim copy of their explanation:

An RNA-Seq dataset from a pilot study of the defense response of
\emph{Arabidopsis} to infection by bacteria. We performed RNA-Seq
experiments on three independent biological samples from each of
the two treatment groups.  The matrix contains the frequencies of
RNA-Seq reads mapped to genes in a reference database. Rows
correspond to genes and columns correspond to independent
biological samples.
}
\details{
This dataset was imported from NBPSeq package and the following
explanation is verbatim copy of their explanation:

We challenged leaves of \emph{Arabidopsis} with the defense-eliciting
\emph{\eqn{\Delta}hrcC} mutant of \emph{Pseudomonas syringae} pathovar
\emph{tomato}
DC3000.  We also infiltrated leaves of \emph{Arabidopsis} with 10mM
MgCl2 as a mock inoculation. RNA was isolated 7 hours after
inoculation, enriched for mRNA and prepared for RNA-Seq. We
sequenced one replicate per channel on the Illumina Genome
Analyzer (http://www.illumina.com).  The length of the RNA-Seq
reads can vary in length depending on user preference and the
sequencing instrument. The dataset used here are derived from a
36-cycle sequencing reaction, that we trimmed to 25mers.  We
used an in-house computational pipeline to process, align, and
assign RNA-Seq reads to genes according to a reference database
we developed for \emph{Arabidopsis}.
}
\usage{data(arab)}
\format{A 26222 by 6 matrix of RNA-Seq read frequencies.}
\examples{
data(arab)
}
\references{
Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
Model for Assessing Differential Gene Expression from RNA-Seq", Statistical
Applications in Genetics and Molecular Biology, 10 (1). 
}
\keyword{datasets}
