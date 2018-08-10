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
ae.input.graphs(gfile = "file.txt", pfile = "file.txt", labelsize = 5)
}
\arguments{
  \item{genotype}{TSV containing SNP information for each strain. First 3 columns are SNP, CHR, BP. Rest of columns are snps for each strain.}
  \item{phenotype}{TSV containing Phenotype information for each individual. Must have column name "Strain" and column name "Value"}
  \item{labelsize}{Size of labels for kinship matrix heatmap. Default is 5.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  ## Set Working Directory
  setwd("~/Auto Emma/sample_files")

  ## Create input graphs
  ae.input.graphs(gfile = "file.txt", pfile = "file.txt", labelsize = 5)
}
