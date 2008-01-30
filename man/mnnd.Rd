\name{mnnd}
\alias{mnnd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mean nearest relative phylogenetic distance }
\description{
  Calculates MNND (mean nearest neighbour or relative phylogenetic distance) for a community 
}
\usage{
mnnd(samp, phy.dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{ Phylocom sample file }
  \item{phy.dist}{ Cophenetic phylogenetic distance matrix }
}
\value{
	Vector of MNND values for each community
}
\references{ Webb et al. 2002. ARES. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\seealso{ \code{\link{nti}} }
