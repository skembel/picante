\name{comdist}
\alias{comdist}
\title{ Calculates inter-community mean pairwise distance }
\description{
  Calculates MPD (mean pairwise distance) separating taxa in two communities, a measure of phylogenetic beta diversity
}
\usage{
comdist(comm, dis, abundance.weighted = FALSE)

}

\arguments{
  \item{comm}{ Community data matrix }
  \item{dis}{ Interspecific distance matrix }
  \item{abundance.weighted}{ Should mean pairwise distances separating species in two communities be weighted by species abundances? (default = FALSE)}
}  
\value{
	Distance object of MPD values separating each pair of communities.
}
\details{
    This function calculates a measure of phylogenetic beta diversity: the expected phylogenetic distance separating two individuals or taxa drawn randomly from different communities.
}
\references{ C.O. Webb, D.D. Ackerly, and S.W. Kembel. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100. }
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{mpd}}, \code{\link{ses.mpd}}  }
\examples{
data(phylocom)
comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)}
\keyword{univar}
