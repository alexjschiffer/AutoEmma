\name{ae.gwl}
\alias{ae.gwl}
\title{Generate Genome Wide Significance Line}
\description{
(This function generates the value of genome-wide significance through a simple bonferroni correction calculated by dividing your original p-value by the number of SNPs, then calculates the -log10 of this number.
}
\usage{
ae.gwl(pval, snps)
}
\arguments{
  \item{pval}{A double}
  \item{snps}{An integer}
}
\value{
The -log10 of \code{pval} divided by \code{snps}
}

\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\examples{
  library(AutoEmma)

  ## Calculate Genome Wide Line
  gwl <- ae.gwl(pval, snps)
}
