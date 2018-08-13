\name{ae.make.map}
\alias{ae.make.map}
\title{Create a Map File Linking Individuals to Strains}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
This function creates a matrix which maps each individual to a particular strain. It both saves a TSV of the map and returns the map as a matrix.
}
\usage{
ae.make.map(file)
}
\arguments{
  \item{file}{A TSV with the first column containing individuals strain name and the first row containing unique strain names.}
}
\value{
  \item{Map}{A binary matrix mapping individuals with strains.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  ## Set Working Directory
  setwd("~/Auto Emma/sample_files")

  map_matrix <- ae.make.map("file_name.txt")
}
