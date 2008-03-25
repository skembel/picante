\name{phylosignal}
\alias{phylosignal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Measure phylogenetic signal }
\description{
  Calculates K statistic of phylogenetic signal as well as P-value based on variance of phylogenetically independent contrasts relative to tip shuffling randomization.
}
\usage{
phylosignal(x, phy, reps = 999, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Trait vector (same order as phy\$tip.label) }
  \item{phy}{ phylo object }
  \item{reps}{ Number of randomizations }
  \item{...}{ Additional arguments passed to pic }
}
\value{
 Data frame with columns...
  \item{K }{K statistic}
  \item{PIC.variance }{Mean observed PIC variance}
  \item{PIC.variance.P }{P-value of observed vs. random variance of PICs}
   \item{PIC.variance.z }{Z-score of observed vs. random variance of PICs}
}
\references{ Blomberg and Garland. 2002. Evolution. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\seealso{ \code{\link{Kcalc}} }
\keyword{univar}
