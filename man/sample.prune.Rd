\name{sample.prune}
\alias{sample.prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prune tree to include species in community }
\description{
  Prune tree to include species in community
}
\usage{
sample.prune(samp, phylo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{ Community data matrix }
  \item{phylo}{ phylo object }
}

\value{
  Returns a phylo object pruned to include only the species present in samp.
}
\author{ Steve Kembel <skembel@berkeley.edu> }
