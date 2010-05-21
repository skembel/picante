\name{pic3}
\alias{pic3}

\title{ Independent contrasts }
\description{
  Modification of pic function from ape to return internal node values as well as contrasts and variance of contrasts.
}
\usage{
pic3(x, phy, scaled = TRUE, var.contrasts = TRUE)
}

\arguments{
  \item{x}{ Trait vector (same order as phy\$tip.label) }
  \item{phy}{ phylo object }
  \item{scaled}{ TRUE = standardize contrasts by branch lengths }
  \item{var.contrasts}{ TRUE = return variance of contrasts }
}
\value{
 Matrix with three columns:
  \item{row.names }{node numbers}
  \item{contrasts }{independent contrast value}
  \item{variance }{variance of contrasts (sum of subtending branch lengths)}
  \item{node.vals }{values calculated at internal nodes by contrast algorithm}
}
\references{ 
Felsenstein J (1985) Phylogenies and the comparative method.  American  Naturalist 125: 1-15.

Garland T, Jr, Harvey PH & Ives AR (1992) Procedures for the analysis of comparative data using phylogenetically independent contrasts.  Systematic  Biology 41: 18-32.
 }
\author{ Original pic function by Emmanuel Paradis <Emmanuel.Paradis@mpl.ird.fr> modified to return node.vals by David Ackerly <dackerly@berkeley.edu> }
\seealso{ \code{\link[ape]{pic}} }
\section{Warning }{ Assumes that trait data are sorted in the same order as phylo\$tip.label }
\examples{
randtree <- rcoal(20)
randtraits <- rTraitCont(randtree)
pic3(randtraits,randtree)
}
\keyword{univar}