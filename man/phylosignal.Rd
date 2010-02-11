\name{phylosignal}
\alias{phylosignal}

\title{ Measure phylogenetic signal }
\description{
  Calculates K statistic of phylogenetic signal as well as P-value based on variance of phylogenetically independent contrasts relative to tip shuffling randomization.
}
\usage{
phylosignal(x, phy, reps = 999, checkdata=TRUE, ...)
}

\arguments{
  \item{x}{ Trait vector (same order as phy\$tip.label) }
  \item{phy}{ phylo object }
  \item{reps}{ Number of randomizations }
  \item{checkdata}{ Check for match between trait and phylogeny taxa labels using \code{\link{match.phylo.data}}? (default=TRUE)}  
  \item{...}{ Additional arguments passed to pic }
}
\value{
 Data frame with columns:
  \item{K }{K statistic}
  \item{PIC.variance }{Mean observed PIC variance}
  \item{PIC.variance.P }{P-value of observed vs. random variance of PICs}
   \item{PIC.variance.z }{Z-score of observed vs. random variance of PICs}
}
\references{ 
Blomberg, S. P., and T. Garland, Jr. 2002. Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology 15:899-910.

Blomberg, S. P., T. Garland, Jr., and A. R. Ives. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57:717-745. }
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{Kcalc}} }
\examples{
randtree <- rcoal(20)
randtraits <- evolve.brownian(randtree)
phylosignal(randtraits[randtree$tip.label],randtree)
}
\keyword{univar}