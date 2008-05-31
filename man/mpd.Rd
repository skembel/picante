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
  \item{dis}{ Interspecific distance matrix }
}

\value{
  Vector of MPD values for each community
}
\references{ Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002. Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\seealso{ \code{\link{ses.mpd}}}
\keyword{univar}
