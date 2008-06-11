\name{pd}
\alias{pd}

\title{ Calculate Faith's Phylogenetic Diversity }
\description{
  Calculate the sum of the total branch lengths for one or multiple samples.
}
\usage{
pd(samp, tree)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{tree}{ A phylo tree object }
}

\value{
	Returns a dataframe of the PD and species richness values for all samples
}

\note{This function does not calculate PD for samples with one or zero species. The data sets need not be species-community data sets but may be any sample
  data set with an associated phylogeny. PD is not statistically independent of species richness, it positively correlates with species richness across samples.
}


\references{ Faith D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10 }
\author{ Matthew Helmus \email{mrhelmus@gmail.com}\cr
         Jonathan Davies \email{davies@nceas.ucsb.edu}
 }
\seealso{ \code{\link{psr}} }

\keyword{univar}
