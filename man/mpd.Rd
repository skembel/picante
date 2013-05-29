\name{mpd}
\alias{mpd}
\title{ Mean pairwise distance }
\description{
  Calculates mean pairwise distance separating taxa in a community
}
\usage{
mpd(samp, dis, abundance.weighted=FALSE)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{dis}{ Interspecific distance matrix }
  \item{abundance.weighted}{ Should mean pairwise distances be weighted by species abundance? (default = FALSE)}  
}

\value{
  Vector of MPD values for each community
}
\references{ Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002. Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505. }
\author{ Steven Kembel <steve.kembel@gmail.com> }
\seealso{ \code{\link{ses.mpd}}}
\examples{
data(phylocom)
mpd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
}
\keyword{univar}
