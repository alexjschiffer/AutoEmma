\name{emma.REML.t}
\alias{emma.REML.t}
\title{Linear mixed model association via t-test with REML estimates.}
\description{
  Performs an efficient linear mixed model association mapping via
  t-test after estimating variance component using REML.
  Progress bar added to AutoEmma version.
}
\usage{
  emma.REML.t (ys, xs, K, Z=diag(ncol(ys)),
  X0 = matrix(1,nrow(ys),1), ngrids=100, llim=-5, ulim=5,
  esp=1e-10, ponly = FALSE, eigen.R0 = NULL, eigen.R1 = NULL)
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
  \item{eig.R0}{Eigenvector from X0, Z and K used in REML estimate. If
    specified, it may avoid redundant computation inside the function}
  \item{eig.R1}{Eigenvector from X1, Z and K used in REML estimate. If
    specified, it may avoid redundant computation inside the function.
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
  \item{ps}{The m by g matrix of p-values between every pair of
    indicator-response variables}
  \item{REMLs}{The g by m matrix of restricted maximum likelihoods}
  \item{stats}{The g by m matrix of t-statistic values}
  \item{dfs}{The m by g matrix of degrees of freedoms}
  \item{vgs}{The m by g matrix of genetic variance components in REML estimates}
  \item{ves}{The m by g matrix of random variance components in REML estimates}

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
    rs <- emma.REML.t(emmadat$ys,emmadat$xs,emmadat$K)

    ## return p-values
    rs
}
