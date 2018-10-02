\name{emma.kinship}
\alias{emma.kinship}
\title{IBS Kinship matrix generation}
\description{
  Generates and returns an IBS kinship matrix from snp data. Progress bar added to AutoEmma version.
}
\usage{
  emma.kinship(xs, method = c("additive","dominant","recessive"), use = c("all","complete.obs"))
}
\arguments{
  \item{xs}{a m by t matrix, where m is number of indicator variables
    (or snps), and n is the number of strains}
  \item{method}{a character string representing the effect of
    heterozygous alleles.}
  \item{use}{a character string giving a method for computing kinship
    coefficient in the presence of missing values. This must be either "all" or "complete.obs"}
}
\details{
  The IBS kinship matrix is computed as simple pairwise genotype similarities.
  Based on the one of additive, dominant, recessive models
  for heterozygous alleles. If "all" option is used, the missing alleles
  are used based on the minor allele frequencies. If "complete.obs"
  option is used, all SNPs containing missing alleles are discarded.
}
\value{
  A t by t matrix containing kinship coefficients between every pair of strains.
}
\author{Hyun Min Kang \email{h3kang@cs.ucsd.edu}}
\references{
  Kang HM, Zaitlen NA, Wade CM, Kirby A, Heckerman D, Daly MJ, and Eskin
  E, Efficient Control of Population Structure in Model Organism
  Association Mapping, Genetics 178:1709-1723, 2008
}
