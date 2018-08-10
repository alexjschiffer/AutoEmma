\name{emma.ML.LRT}
\alias{emma.ML.LRT}
\title{Linear mixed model association via Likelihood Ratio test with ML estimates.}
\description{
  Performs an efficient linear mixed model association mapping via
  Likelihood Ratio test (LR) after estimating variance component using REML.
  Progress bar added to AutoEmma version.
}
\usage{
  emma.REML.LRT (ys, xs, K, Z=diag(ncol(xs)),
  X0 = matrix(1,nrow(xs),1), ngrid=100, llim=-5, ulim=5,
  esp=1e-10, ponly = FALSE)
}
\arguments{
  \item{ys}{A g by n matrix, where g is the number of response variables
    (or phenotypes), and n is the number of individuals}
  \item{xs}{A m by t matrix, where m is number of indicator variables
    (or snps), and n is the number of strains}
  \item{K}{A t by t matrix of kinship coefficients, representing the
    pairwise genetic relatedness between strains}
  \item{Z}{A n by t incidence matrix mapping each individual to a
    strain. If this is NULL, n and t should be equal and an identity
    matrix replace Z}
  \item{X0}{A n by p matrix of fixed effects variables, where p is the
    number of fixed effects including mean and other confounding variables}
  \item{ngrids}{Number of grids to search optimal variance component}
  \item{llim}{Lower bound of log ratio of two variance components}
  \item{ulim}{Upper bound of log ratio of two variance components}
  \item{esp}{Tolerance of numerical precision error}
  \item{ponly}{Returns p-value matrix only if TRUE}
  \item{eig.L}{Eigenvector from Z and K used in ML estimate. If
    specified, it may avoid redundant computation inside the function}
  \item{eig.R0}{Eigenvector from X0, Z and K used in REML estimate. If
    specified, it may avoid redundant computation inside the function}
  \item{eig.R1}{Eigenvector from X1, Z and K used in REML estimate. If
    specified, it may avoid redundant computation inside the function
    Valid only when m=1}
}
\details{
  The following criteria must hold; otherwise an error occurs
  - [# cols in ys] == [# rows in Z] == [# rows in X0]
  - [# cols in xs] == [# cols in Z] == [# rows in K] == [# cols in K]
  - rowSums(Z) should be a vector of ones
  - colSums(Z) should not contain zero elements
  - K must be a positive semidefinite matrix
}
\value{
  A list containing:
  \item{ps}{The g by m matrix of p-values between every pair of
    indicator-response variables}
  \item{stats}{The g by m matrix of chi-square statistic values}
  \item{ML0s}{The g by m matrix of log maximum likelihood under null hypothesis}
  \item{ML1s}{The g by m matrix of log maximum likliehood under
    alternative hypothesis}
  \item{vgs}{The g by m matrix of genetic variance components under null
    hypothesis}
  \item{ves}{The g by m matrix of random variance components under null
    hypothesis}

  if ponly is TRUE, only ps is return as matrix form.
}
\references{
  Kang HM, Zaitlen NA, Wade CM, Kirby A, Heckerman D, Daly MJ, and Eskin
  E, Efficient Control of Population Structure in Model Organism
  Association Mapping, Genetics 178:1709-1723, 2008
}
\author{Hyun Min Kang \email{h3kang@cs.ucsd.edu}}
\examples{
    ## Load data
    data(emmadat)

    ## Run EMMA
    rs <- emma.ML.LRT(emmadat$ys,emmadat$xs,emmadat$K)

    ## return p-values
    rs$ps
}
