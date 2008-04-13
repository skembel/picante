\name{mpd}
\alias{mpd}
\title{ Mean pairwise distance }
\description{
  Calculates mean pairwise distance for a community
}
\usage{
mpd(samp, dis)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{dis}{ Interspecific (e.g. phylogenetic) distance matrix }
}

\value{
  Vector of MPD values for each community
}
\references{ Webb et al. 2002. ARES. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\seealso{ \code{\link{ses.mpd}}}
\keyword{univar}
