\name{ae.input.graphs}
\alias{ae.input.graphs}
\title{Create Input File Diagnostic Plots}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
Creates kinship matrix heatmap and population structure cladogram for genotype diagnostics. Creates phenotype distribution and average phenotype by strain for phenotype diagnostics.
}
\usage{
ae.input.graphs(gfile = "file.txt", pfile = "file.txt", stdev = FALSE, labelsize = 10)
}
\arguments{
  \item{gfile}{TSV or CSV containing SNP information for each strain. First 3 columns are SNP, CHR, BP. Rest of columns are snps for each strain.}
  \item{pfile}{TSV or CSV containing Phenotype information for each individual. Must have column name "Strain" and column name "Value"}
  \item{stdev}{Boolean indicated whether to use stdev or sem for error bars. Default is FALSE, meaning SEM will be used.}
  \item{labelsize}{Size of labels for kinship matrix heatmap. Default is 10.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  ## Set Working Directory
  setwd("~/Auto Emma/sample_files")

  ## Create input graphs
  ae.input.graphs(gfile = "file.txt", pfile = "file.txt", stdev = FALSE, labelsize = 10)
}
