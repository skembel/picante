\name{expected.pd}
\alias{expected.pd}
\alias{variance.pd}
\alias{ead}

\title{ Expected PD, PD Variance, and Edge Abundance Distribution of a phylogeny }

\description{ Calculates the expected phylogenetic diversity (Faith's PD) and variance of PD under binomial sampling with a fixed probability of each tip being sampled, and the Edge-length Abundance Distribution of a phylogeny. }

\usage{
expected.pd(phy)
variance.pd(phy, upper.bound=TRUE)
ead(phy)
}

\arguments{
  \item{phy}{ phylo object }
  \item{upper.bound}{ Calculate upper bound of PD variance? (default = TRUE) }
}

\value{
  \item{n}{ Expected Number of tips sampled }
  \item{expected.pd}{ Expected PD for a given n }
  \item{variance.pd}{ Variance of PD for a given n }
  \item{num.children}{ Number of tips descended from an edge }
  \item{edge.length}{ Total phylogenetic edge length for a given number of tips descended from an edge }      
}

\details{The function \code{expected.pd} calculates the expected phylogenetic diversity (Faith's PD - total branch length) for all subsets of a phylogeny, based on an analytic solution for expected PD.

The function \code{variance.pd} additionally calculates the variance of expected PD for all subsets of a phylogeny, based on an analytic solution for expected PD. If argument upper.bound=TRUE, a fast solution for the upper bound of the variance is returned. Otherwise, the exact solution for the variance is returned. Note that the exact solution is much slower than the upper bound solution.

The function \code{ead} calculates the edge abundance distribution (EAD), the length of edges with different numbers of descendant tips.
}

\references{ J.P. O'Dwyer, S.W. Kembel, and J.L. Green. 2012. Phylogenetic Diversity Theory Sheds Light on the Structure of Microbial Communities. PLoS Comput Biol 8(12): e1002832. }

\author{ Steven Kembel <steve.kembel@gmail.com> and James O'Dwyer <jodwyer@santafe.edu> }

\seealso{ \code{\link{pd}} }

\examples{
randtree <- rcoal(300)
randtree.pd.ub <- variance.pd(randtree, upper.bound=TRUE)
randtree.pd.exact <- variance.pd(randtree, upper.bound=FALSE)
plot(expected.pd(randtree), xlab="Number of tips",
    ylab="Phylogenetic diversity (PD)", type="l", log="xy")
lines(randtree.pd.exact$expected.pd+1.96*sqrt(randtree.pd.exact$variance.pd), lty=2)
lines(randtree.pd.exact$expected.pd-1.96*sqrt(randtree.pd.exact$variance.pd), lty=2)
lines(randtree.pd.ub$expected.pd+1.96*sqrt(randtree.pd.ub$variance.pd), lty=3)
lines(randtree.pd.ub$expected.pd-1.96*sqrt(randtree.pd.ub$variance.pd), lty=3)
legend("bottomright", lty=c(1,2,3), legend=c("Expected PD",
    "95 percent CI (exact)","95 percent CI (upper bound)"))
}
\keyword{univar}
