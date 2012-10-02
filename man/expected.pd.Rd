\name{expected.pd}
\alias{expected.pd}
\alias{ead}

\title{ Expected PD and Edge Abundance Distribution of a phylogeny }

\description{ Calculates the expected phylogenetic diversity (Faith's PD) under binomial sampling with a fixed probability of each tip being sampled, and the Edge-length Abundance Distribution of a phylogeny. }

\usage{
expected.pd(phy)
ead(phy)
}

\arguments{
  \item{phy}{ phylo object }
}

\value{
  \item{n}{ Expected Number of tips sampled }
  \item{expected.pd}{ Expected PD for a given n }
  \item{num.children}{ Number of tips descended from an edge }
  \item{edge.length}{ Total phylogenetic edge length for a given number of tips descended from an edge }      
}

\details{The function \code{expected.pd} calculates the expected phylogenetic diversity (Faith's PD - total branch length) for all subsets of a phylogeny, based on an analytic solution for expected PD.

The function \code{ead} calculates the edge abundance distribution (EAD), the length of edges with different numbers of descendant tips.
}

\references{ O'Dwyer, Kembel, and Green. 2012. Phylogenetic diversity theory sheds light on the structure of microbial communities. }

\author{ Steven Kembel <skembel@uoregon.edu> and James O'Dwyer <jodwyer@santafe.edu> }

\seealso{ \code{\link{pd}} }

\examples{
randtree <- rcoal(100)
plot(expected.pd(randtree), xlab="Number of tips", ylab="Expected PD", type="l")
}
\keyword{univar}
