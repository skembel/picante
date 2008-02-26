\name{species.dist}
\alias{species.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute various indices of species co-occurrence }
\description{
  Compute various indices of species co-occurrence
}
\usage{
species.dist(samp, metric = c("cij", "jaccard", "roij"))
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{metric}{ Distance metric to use (COij, Jaccard, ROij) }
}

\value{
  dist object
}

\author{ Steve Kembel <skembel@berkeley.edu> }

\seealso{ \code{\link{cij}}, \code{\link{roij}},\code{\link[vegan]{vegdist}} }
\keyword{univar}
