\name{ae.manhattan}
\alias{ae.manhattan}
\title{Creates and Save a Manhattan Plot}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
This function creates and save a manhattan plot.
}
\usage{
  ae.manhattan(df, file_name, point_size = 1.2, colors = c("royalblue3", "gray18"), title = "Manhattan Plot", sl = 4.5, gwl = 7, exclude = 0.5)
}
\arguments{
  \item{df}{A dataframe with at least "SNP", "CHR", "BP", "P".}
  \item{file_name}{A string to name the resulting png file.}
  \item{point_size}{An integer describing the size of each point.}
  \item{colors}{A vector describing the colors for the plot. Default: c("royalblue3", "gray18")}
  \item{sl}{A double describing the position of the suggestive line. Default: 4.5}
  \item{gwl}{A double describing the position of the genome-wide significance line. Default: 7}
  \item{title}{A string to name the chromosome plots. Default: "Manhattan Plot"}
  \item{exclude}{A double describing the cutoff -log10(p) to exclude points below. Speeds up graphing.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  ## Set Working Directory
  setwd("~/Auto Emma/sample_files")

  ## Calculate Genome Wide Line
  line <- ae.gwl(pval, snps)

  ## Make Manhattan Plot
  ae.manhattan(df, file_name, point_size = 1.2, colors = c("royalblue3", "gray18"), title = "Manhattan Plot", sl = 4.5, gwl = 7, exclude = 0.5)
}
