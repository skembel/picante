\name{comdistnt}
\alias{comdistnt}
\alias{comdistnn}
\title{ Calculates inter-community mean nearest taxon distance }
\description{
  Calculates MNTD (mean nearest taxon distance) separating taxa in two communities, a measure of phylogenetic beta diversity
}
\usage{
comdistnt(comm, dis, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
}

\arguments{
  \item{comm}{ Community data matrix }
  \item{dis}{ Interspecific distance matrix }
  \item{abundance.weighted}{ Should mean nearest taxon distances from each species to species in the other community be weighted by species abundance? (default = FALSE)}
  \item{exclude.conspecifics}{ Should conspecific taxa in different communities be exclude from MNTD calculations? (default = FALSE) }  
}
\value{
	Distance object of MNTD values separating each pair of communities.
}
\details{
    This metric has also been referred to as MNND (mean nearest neighbour distance).
    
    This function calculates a measure of phylogenetic beta diversity: the  average phylogenetic distance to the most similar taxon or individual in the other community for taxa or individuals in two communities.
}
\references{ C.O. Webb, D.D. Ackerly, and S.W. Kembel. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100. }
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{mntd}}, \code{\link{ses.mntd}} }
\examples{
data(phylocom)
comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)}
\keyword{univar}
