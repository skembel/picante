\name{matchSpeciesMatrix}
\alias{matchSpeciesMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ensure rows and columns of two community data matrices match }
\description{
  Ensure rows and columns of two community data matrices match
}
\usage{
matchSpeciesMatrix(x, y)
}

\arguments{
  \item{x}{ Community data matrix 1}
  \item{y}{ Community data matrix 2}
}
\value{
  Community data matrix 2 but with rows and columns for all sites and species present in community data matrix 1 present (even if those species/sites not present in matrix 2)
}
\author{ Steve Kembel <skembel@berkeley.edu> }
