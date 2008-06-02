\name{utility}
\alias{df2vec}
\alias{internal2tips}
\alias{node.age}
\alias{pic.variance}
\alias{sortColumns}
\alias{sortRows}
\alias{taxaShuffle}
\title{ Picante utility functions }
\description{
  Picante utility functions for tree and data manipulation
}
\usage{
df2vec(x, colID=1)
internal2tips(phy, int.node, return.names = FALSE)
node.age(phy)
pic.variance(x, phy, scaled = TRUE)
sortColumns(x)
sortRows(x)
taxaShuffle(x)
}

\arguments{
    \item{phy}{ phylo object }
    \item{x}{ A data.frame, matrix or dist object }
    \item{colID}{ Numeric or character ID of column to include }
    \item{int.node}{ internal node number }
    \item{return.names}{ TRUE or FALSE }
    \item{scaled}{ Scale contrasts by branch length }
    \item{...}{Additional arguments}
}

\value{
  \item{df2vec}{A named vector}
  \item{internal2tips}{Vector of tips descended from a node}
  \item{node.age}{Phylo object with phylo\$ages vector of node ages corresponding to phylo\$edge}
  \item{pic.variance}{Variance of independent contrasts}
  \item{sortColumns}{A data.frame or matrix with columns sorted by name}
  \item{sortRows}{A data.frame or matrix with rows sorted by name}
    \item{taxaShuffle}{Matrix with taxa names shuffled}
}

\details{
  Various utility functions for manipulating trees, data, etc.
}

\author{ Steven Kembel <skembel@berkeley.edu>, Peter Cowan <pdc@berkeley.edu>, David Ackerly <dackerly@berkeley.edu> }
\keyword{manip}
