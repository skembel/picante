\name{mntd}
\alias{mntd}
\alias{mnnd}
\title{ Mean nearest taxon distance }
\description{
  Calculates MNTD (mean nearest taxon distance) for taxa in a community
}
\usage{
mntd(samp, dis, abundance.weighted=FALSE)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{dis}{ Interspecific distance matrix }
  \item{abundance.weighted}{ Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)}
}
\value{
	Vector of MNTD values for each community.
}
\details{
    This metric has also been referred to as MNND (mean nearest neighbour distance), and the function was named \code{mnnd} in picante versions < 0.7.
}
\references{ Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002. Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505. }
\author{ Steven Kembel <steve.kembel@gmail.com> }
\seealso{ \code{\link{ses.mntd}} }
\examples{
data(phylocom)
mntd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
}
\keyword{univar}
