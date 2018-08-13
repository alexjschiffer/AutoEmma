\name{ae.make.regional}
\alias{ae.make.regional}
\title{Make Regional Chromosome Plots}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
  Creates regional manhattan plots for each chromosome and saves them as a single PDF.
}
\usage{
  ae.make.regional(resultsdata = dataframe, filename = "filename.pdf"), chromosomes = c(1:20), sl = 4.5, gwl = 7, title = "Regional Manhattan Plot")
}
\arguments{
  \item{resultsdata}{A dataframe with at least "SNP", "CHR", "BP", "P".}
  \item{filename}{A string to name the resulting pdf file.}
  \item{chromosomes}{A vector describing the number of chromosomes. Default: c(1:20)}
  \item{sl}{A double describing the position of the suggestive line. Default: 4.5}
  \item{gwl}{A double describing the position of the genome-wide significance line. Default: 7}
  \item{title}{A string to name the chromosome plots. Default: "Regional Manhattan Plot"}
}
\value{
  \item{None}{Graphs are saved to folder. No value is returned.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  # Import results data
  dataframe <- ae.import("resultsdata.csv")

  # Make regional plots
  ae.make.regional(resultsdata = dataframe, filename = "filename.pdf"), chromosomes = c(1:20), sl = 4.5, gwl = 7, title = "Regional Manhattan Plot")
}
