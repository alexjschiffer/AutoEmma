\name{ae.qq}
\alias{ae.qq}
\title{Creates and Save a Q-Q Plot}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
This function creates and save a Q-Q plot.
}
\usage{
  ae.qq(df, file_name, title = "Q-Q Plot")
}
\arguments{
  \item{df}{A dataframe with at least "SNP", "CHR", "BP", "P".}
  \item{file_name}{A string to name the resulting png file.}
  \item{title}{A string to name the chromosome plots. Default: "Manhattan Plot"}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  ## Set Working Directory
  setwd("~/Auto Emma/sample_files")

  ## Make QQ Plot
  ae.qq(df, file_name, title = "Q-Q Plot")
}
