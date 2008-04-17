\name{Kcalc}
\alias{Kcalc}

\title{ K statistic of phylogenetic signal }
\description{
  Calculates K statistic of phylogenetic signal
}
\usage{
Kcalc(x, phy)
}

\arguments{
  \item{x}{ Vector or data.frame of trait data (in phylo\$tip.label order) }
  \item{phy}{ phylo object }
}

\value{
  \item{K }{K statistic}
}
\references{ Blomberg, S. P., T. Garland, Jr., and A. R. Ives. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57:717-745.}
\author{ Simon Blomberg <s.blomberg1@uq.edu.au> and David Ackerly <dackerly@berkeley.edu> }

\seealso{ \code{\link{phylosignal}} }
\keyword{univar}
