\name{match.phylo.data}
\alias{match.phylo.data}
\alias{match.phylo.comm}

\title{ Match taxa in phylogeny and data }
\description{
  These functions compare taxa present in phylogenies with community or trait data sets, pruning and sorting the two kinds of data to match one another for subsequent analysis.
}
\usage{
match.phylo.comm(phy, comm)
match.phylo.data(phy, data)
}

\arguments{
  \item{ phy }{ A phylogeny object of class phylo }
  \item{ comm }{ Community data matrix }
  \item{ data }{ A data object - a vector (with names matching phy) or a data.frame or matrix (with row names matching phy) }
}

\value{
  A list containing the following elements, pruned and sorted to match one another:
  \item{phy}{ A phylogeny object of class phylo }
  \item{comm}{ Community data matrix }
  \item{data}{ A data object (vector, data.frame or matrix) }  
}
\details{
A common pitfall in comparative analyses in R is that taxa labels are assumed to match between phylogenetic and other data sets. These functions prune a phylogeny and community or trait data set to match one another, reporting taxa that are missing from one data set or the other.

Taxa names for phylogeny objects are taken from the phylogeny's tip labels. Taxa names for community data are taken from the column names. Taxa names for trait data are taken from the element names (vector) or row names (data.frame or matrix).

If community data lack taxa names, the function will issue a warning and no result will be returned, since the community-phylogenetic analyses in \code{picante} require named taxa in the community data set.

If trait data lack names, a warning is issued and the data are assumed to be sorted in the same order as the phylogeny's tip labels.

These utility functions are used by several of the functions that assume taxa labels in phylogeny and data match, including \code{\link{Kcalc}}, \code{\link{phylosignal}}, and \code{\link{raoD}}.
}
\author{ Steven Kembel <skembel@uoregon.edu> }
\examples{
data(phylocom)
match.phylo.comm(phylocom$phylo, phylocom$sample)
match.phylo.data(phylocom$phylo, phylocom$traits[1:10,])
}
\keyword{univar}
