\name{ae.make.manhattan}
\alias{ae.make.manhattan}
\title{Creates and Save a Manhattan Plot}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
This function creates and save a manhattan plot.
}
\usage{
  ae.manhattan <- function(resultsdata = dataframe, filename = "manhattan.png", colors = c("royalblue3", "gray18"), sl = 4.5, gwl = 7, title = "Manhattan Plot")
}
\arguments{
  \item{resultsdata}{A dataframe with at least "SNP", "CHR", "BP", "P".}
  \item{filename}{A string to name the resulting png file.}
  \item{colors}{A vector describing the colors for the plot. Default: c("royalblue3", "gray18")}
  \item{sl}{A double describing the position of the suggestive line. Default: 4.5}
  \item{gwl}{A double describing the position of the genome-wide significance line. Default: 7}
  \item{title}{A string to name the chromosome plots. Default: "Manhattan Plot"}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  ## Set Working Directory
  setwd("~/Auto Emma/sample_files")

  ## Calculate Genome Wide Line
  line <- ae.gwl(pval, snps)

  ## Make Manhattan Plot
  ae.manhattan <- function(resultsdata = dataframe, filename = "manhattan.png", colors = c("royalblue3", "gray18"), sl = 4.5, gwl = line, title = "Manhattan Plot")
}
