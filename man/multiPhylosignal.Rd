\name{multiPhylosignal}
\alias{multiPhylosignal}
\title{ Calculates phylogenetic signal for data.frame of traits }
\description{
  Calculates phylogenetic signal for data.frame of traits. Traits may have missing values in which case the tree will be pruned prior to calculating phylogenetic signal for each trait.
}
\usage{
multiPhylosignal(x, phy, ...)
}
\arguments{
  \item{x}{ Data frame of trait data (traits in columns) with row names corresponding to tip.labels}
  \item{phy}{ phylo object }
  \item{...}{ Additional arguments to \link{phylosignal} }
}
\value{
  Returns a data frame with phylogenetic signal results for each trait
}
\section{Warning }{ Assumes that trait data are sorted in the same order as phylo tip.label }
\author{ Steve Kembel <skembel@berkeley.edu> }
\keyword{univar}
