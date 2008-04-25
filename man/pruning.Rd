\name{prune.sample}
\alias{prune.sample}
\alias{prune.missing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prune tree to match community data or trait data }
\description{
  Prune a phylogenetic tree to include only species present in a community data set or with trait data}
\usage{
prune.sample(samp, phylo)
prune.missing(x, phylo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{ phylo object }
  \item{samp}{ Community data matrix }
  \item{x}{ trait data }
}

\value{
  Returns a pruned phylo object
}
\author{ Steve Kembel <skembel@berkeley.edu> }
\keyword{manip}
