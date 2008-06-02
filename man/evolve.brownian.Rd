\name{evolve.brownian}
\alias{evolve.brownian}
\title{ Trait evolution under Brownian motion }
\description{
  Trait evolution under Brownian motion model
}
\usage{
evolve.brownian(phy, value = 0, var = 1)
}

\arguments{
  \item{phy}{ phylo object }
  \item{value}{ mean }
  \item{var}{ variance }
}
\value{
  Vector of trait values with names corresponding to phylo\$tip.label
}
\author{ David Ackerly <dackerly@berkeley.edu> and Steve Kembel <skembel@berkeley.edu> }
\keyword{datagen}
