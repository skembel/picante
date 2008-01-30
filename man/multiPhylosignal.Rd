\name{multiPhylosignal}
\alias{multiPhylosignal}
\title{ Calculates phylogenetic signal for data.frame of traits }
\description{
  Calculates phylogenetic signal for data.frame of traits. Traits may have missing values in which case the tree will be pruned prior to calculating phylogenetic signal for each trait.
}
\usage{
multiPhylosignal(x,tree,...)
}
\arguments{
  \item{x}{ Data frame of trait data with row names corresponding to tip.labels }
  \item{tree}{ phylo object }
  \item{...}{ Additional arguments to phylosignal }
}
\value{
  Returns a data frame with phylogenetic signal results for each trait
}
\author{ Steve Kembel <skembel@berkeley.edu> }
