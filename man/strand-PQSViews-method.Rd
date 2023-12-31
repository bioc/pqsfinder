% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PQSViews-class.R
\name{strand,PQSViews-method}
\alias{strand,PQSViews-method}
\title{Get PQS strand vector}
\usage{
\S4method{strand}{PQSViews}(x)
}
\arguments{
\item{x}{PQSViews object.}
}
\value{
Strand vector.
}
\description{
Get PQS strand vector
}
\examples{
pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
strand(pqs)

}
