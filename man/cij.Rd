\name{cij}
\alias{cij}
\title{ Schoener's pairwise index of species co-occurrence }
\description{
  Calculates Schoener's pairwise index of species co-occurrence for a community matrix
}
\usage{
cij(x)
}
\arguments{
  \item{x}{ Community data matrix }
}
\value{
  Return an object of class dist containing Schoener's index of co-occurrence for all species pairs in x.
}
\references{ Schoener }
\author{ Steve Kembel <skembel@berkeley.edu> }
\keyword{univar}
