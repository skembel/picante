\name{Kcalc}
\alias{Kcalc}

\title{ K statistic of phylogenetic signal }
\description{
  Calculates K statistic of phylogenetic signal
}
\usage{
Kcalc(x, phy, checkdata=TRUE)
}

\arguments{
  \item{x}{ Vector or data.frame of trait data (in phylo\$tip.label order) }
  \item{phy}{ phylo object }
  \item{checkdata}{ Check for match between trait and phylogeny taxa labels using \code{\link{match.phylo.data}}? (default=TRUE)}
}

\value{
  \item{K }{K statistic}
}
\references{ Blomberg, S. P., and T. Garland, Jr. 2002. Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology 15:899-910.

Blomberg, S. P., T. Garland, Jr., and A. R. Ives. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57:717-745.}

\author{ Simon Blomberg <s.blomberg1@uq.edu.au> and David Ackerly <dackerly@berkeley.edu> }

\seealso{ \code{\link{phylosignal}} }
\examples{
randtree <- rcoal(20)
randtraits <- rTraitCont(randtree)
Kcalc(randtraits[randtree$tip.label],randtree)
}
\keyword{univar}
