\name{ae.run}
\alias{ae.run}
\title{Automatic Efficient Mixed Model Association}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
  Imports provided files, runs EMMA (UCLA) with either REML or ML, then saves figures and results data to a results folder.
}
\usage{
ae.run(genotype = "", phenotype = "", kinship = NULL, map = NULL, method = "REML", use_snps = "all")
}
\arguments{
  \item{genotype}{TSV containing SNP information for each strain. First 3 columns are SNP, CHR, BP. Rest of columns are snps for each strain.}
  \item{phenotype}{TSV containing Phenotype information for each individual. Must have column name "Strain" and column name "Value"}
  \item{kinship}{A matrix indicating strain x strain relatedness. Default: NULL (function will generate one)}
  \item{map}{Matrix of ones and zeros matching each individual to a strain. Default: NULL}
  \item{method}{A string, either "REML" or "ML". Default: "REML"}
  \item{gwl}{A double. Default: 6}
  \item{use_snps}{a string, either "all" or "complete.obs". Default: "all"}
}
\value{
  \item{results}{A data table containing snp information and p-values for those snps.}
}
\examples{
  # Load AutoEmma R Package
  library(AutoEmma)

  # Set Working Directory
  setwd("~/Auto Emma/sample_files")

  # Create genome-wide line and map file
  line <- ae.gwl(0.05, 108500)
  map1 <- ae.make.map(pheno.txt)

  # Run AutoEMMA
  results <- ae.run(genotype = "", phenotype = "", kinship = NULL, map = NULL, method = "REML", use_snps = "all")
}
