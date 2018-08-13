\name{ae.import}
\alias{ae.import}
\title{Import CSV and TSV Files}
\author{
  Alex Schiffer \email{alex.schiffer@childrens.harvard.edu}
}
\description{
  Imports files using the readr package (.txt and .tab imported as TSV; .csv imported as CSV). Provides an option for loading only unique rows, which can be helpful for removing duplicate entries.
}
\usage{
  ae.import(input = "filename.txt", unique = FALSE)
}
\arguments{
  \item{file}{TSV or CSV file with headers.}
  \item{unique}{Boolean indicating whether to load all rows or only unique rows.}
}
\value{
  \item{results}{A data table with column names.}
}
\examples{
  ae.import(input = "filename.txt", unique = FALSE)
}
