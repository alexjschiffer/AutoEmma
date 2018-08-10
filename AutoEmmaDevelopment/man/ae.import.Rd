\name{ae.import}
\alias{ae.import}
\title{Import}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
  Imports provided files (.txt and .tab imported as TSV; .csv imported as CSV).
}
\usage{
  ae.import("filename.txt")
}
\arguments{
  \item{file}{TSV or CSV file with headers}
}
\value{
  \item{results}{A data table with column names.}
}
\examples{
  ae.import("filename.csv")
}
