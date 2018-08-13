\name{ae.recode}
\alias{ae.recode}
\title{Recodes a SNP file from ATCG to binary}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
This function converts a file from ATCGN to binary. Lowercase are converted to Uppercase first. N are converted to NA.
}
\usage{
ae.recode(dataframe)
}
\arguments{
  \item{dataframe}{A data frame containing snps in atcg format.}
}
\value{
  \item{recoded}{A data table in the recoded format. A TSV is also saved to the working directory.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  # Import results data
  recoded <- ae.recode(dataframe)
}
