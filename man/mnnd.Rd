\name{mnnd}
\alias{mnnd}
\title{ Mean nearest neighbour distance }
\description{
  Calculates MNND (mean nearest neighbour distance) for taxa in a community
}
\usage{
mnnd(samp, dis)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{dis}{ Interspecific distance matrix }
}
\value{
	Vector of MNND values for each community.
}
\references{ Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002. Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505. }
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{ses.mnnd}} }
\keyword{univar}
