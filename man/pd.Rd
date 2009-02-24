\name{pd}
\alias{pd}

\title{ Calculate Faith's Phylogenetic Diversity }
\description{
  Calculate the sum of the total phylogenetic branch length for one or multiple samples.
}
\usage{
pd(samp, tree, include.root=TRUE)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{tree}{ A phylo tree object }
  \item{include.root}{ Should the root node be included in all PD calculations (default = TRUE)}
}

\value{
	Returns a dataframe of the PD and species richness (SR) values for all samples
}

\note{
The data sets need not be species-community data sets but may be any sample data set with an associated phylogeny. PD is not statistically independent of species richness, it positively correlates with species richness across samples.

If the root is to be included in all calculations of PD (\code{include.root=TRUE}), the tree must be rooted and ultrametric. Single-species samples will be assigned a PD value equal to the distance from the root to the present.

If the root is not included in all calculations by default (\code{include.root=FALSE}), the tree need not be ultrametric or rooted, but in the case of single-species samples the PD will be equal to NA and a warning will be issued.

}
\section{Warning }{
If the root is to be included in all calculations (\code{include.root=TRUE}), the PD of all samples will include the branch length connecting taxa in those samples and the root of the supplied tree. The root of the supplied tree may not be spanned by any taxa in the sample. If you want the root of your tree to correspond to the most recent ancestor of the taxa actually present in your sample, you should prune the tree before running \code{pd}:

\code{prunedTree <- prune.sample(sample,tree)}
}
\references{ Faith D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10. }
\author{ Matthew Helmus \email{mrhelmus@gmail.com}, Jonathan Davies \email{davies@nceas.ucsb.edu}, Steven Kembel \email{skembel@uoregon.edu}}
\seealso{ \code{\link{psr}} }

\keyword{univar}
