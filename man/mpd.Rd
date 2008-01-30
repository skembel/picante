\name{mpd}
\alias{mpd}
\title{ Mean pairwise phylogenetic distance }
\description{
  Calculates mean pairwise phylogenetic distance for a community
}
\usage{
mpd(samp, phy.dist)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{phy.dist}{ Cophenetic phylogenetic distance matrix }
}

\value{
  Vector of MPD values for each community
}
\references{ Webb et al. 2002. ARES. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\seealso{ \code{\link{nri}}}
