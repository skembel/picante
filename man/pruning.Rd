\name{prune.sample}
\alias{prune.sample}
\alias{prune.missing}

\title{ Prune tree to match community data or trait data }
\description{
  Prune a phylogenetic tree to include only species present in a community data set or with non-missing trait data}
\usage{
prune.sample(samp, phylo)
prune.missing(x, phylo)
}

\arguments{
  \item{phylo}{ phylo object }
  \item{samp}{ Community data matrix }
  \item{x}{ Vector of trait data }
}

\value{
  Returns a pruned phylo object
}
\author{ Steven Kembel <skembel@uoregon.edu> }
\keyword{manip}
