\name{pruneMissing}
\alias{pruneMissing}
\title{ Prune tips with missing data from a tree }
\description{
  Prune tips with missing data from a tree
}
\usage{
pruneMissing(x,tree)
}
\arguments{
  \item{x}{ Vector of trait data with names corresponding to tip.labels }
  \item{tree}{ phylo object }
}
\value{
  Returns a list
  \item{tree}{ Pruned phylo object }
  \item{data}{ Corresponding trait data sorted in tip.label order }  
}
\author{ Steve Kembel <skembel@berkeley.edu> }
