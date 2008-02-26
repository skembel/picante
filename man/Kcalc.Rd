\name{Kcalc}
\alias{Kcalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ K statistic of phylogenetic signal }
\description{
  Calculates K statistic of phylogenetic signal
}
\usage{
Kcalc(x, phy)
}

\arguments{
  \item{x}{ Vector of trait data (in phylo\$tip.label order) }
  \item{phy}{ phylo object }
}

\value{
  \item{K }{K statistic}
}
\references{ Blomberg and Garland 2002 Evolution}
\author{ Simon Blomberg <s.blomberg1@uq.edu.au> and David Ackerly <dackerly@berkeley.edu> }

\seealso{ \code{\link{phylosignal}} }
