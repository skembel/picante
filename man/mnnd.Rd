\name{mnnd}
\alias{mnnd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mean nearest neighbour distance }
\description{
  Calculates MNND (mean nearest neighbour distance) for a community
}
\usage{
mnnd(samp, dis)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{ Phylocom sample file }
  \item{dis}{ Interspecific (e.g. phylogenetic) distance matrix }
}
\value{
	Vector of MNND values for each community
}
\references{ Webb et al. 2002. ARES. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\seealso{ \code{\link{ses.mnnd}} }
\keyword{univar}
