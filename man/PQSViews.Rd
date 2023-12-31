% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PQSViews-class.R
\name{PQSViews}
\alias{PQSViews}
\title{PQSViews class constructor}
\usage{
PQSViews(
  subject,
  start,
  width,
  strand,
  score,
  density,
  max_scores,
  nt,
  nb,
  nm,
  rl1,
  rl2,
  rl3,
  ll1,
  ll2,
  ll3
)
}
\arguments{
\item{subject}{DNAString or RNAString object.}

\item{start}{Start positions.}

\item{width}{Lengths.}

\item{strand}{Strand specifications.}

\item{score}{Scores.}

\item{density}{Numbers of PQS overlapping at each position in \code{subject}.}

\item{max_scores}{Score of the best PQS found at each position.}

\item{nt}{Tetrad numbers.}

\item{nb}{Bulge counts.}

\item{nm}{Mismatch counts.}

\item{rl1}{Run 1 lengths.}

\item{rl2}{Run 2 lengths.}

\item{rl3}{Run 3 lengths.}

\item{ll1}{Loop 1 lengths.}

\item{ll2}{Loop 2 lengths.}

\item{ll3}{Loop 3 lengths.}
}
\value{
PQSViews object.
}
\description{
User friendly constructor for PQSViews class representing potential
quadruplex forming sequences (PQS). PQSViews is a subclass of
\code{\link{XStringViews}} class and adds two more slots to store
PQS density and PQS score distribution.
}
\details{
Use \code{\link{elementMetadata}} function to get extra PQS features
like number of tetrads, bulges, mismatches or loop lengths.
}
\examples{
pv <- PQSViews(DNAString("GGTGGTGGTGG"), 1, 11, "+", 33, as.integer(rep(1, 11)),
               as.integer(rep(33, 11)), 2, 0, 0, 2, 2, 2, 1, 1, 1)
start(pv)
width(pv)
strand(pv)
score(pv)
density(pv)
maxScores(pv)
elementMetadata(pv)

}
