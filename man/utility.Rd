\name{utility}
\alias{t2p2t}
\alias{df2vec}
\alias{internal2tips}
\alias{node.age}
\alias{phylo2phylog}
\alias{sortColumns}
\alias{sortRows}
\title{ Picante utility functions }
\description{
  Picante utility functions for tree and data manipulation
}
\usage{
df2vec(x, colID=1)
internal2tips(phy, int.node, return.names = FALSE)
node.age(phy)
phylo2phylog(phy, ...)
sortColumns(x)
sortRows(x)
t2p2t(phy, species)
}

\arguments{
    \item{phy}{ phylo object }
    \item{species}{ Vector of species names }
    \item{x}{ A data.frame or matrix }
    \item{colID}{ Numeric or character ID of column to include }
  \item{int.node}{ internal node number }
  \item{return.names}{ TRUE or FALSE }
  \item{...}{Additional arguments}
}

\value{
  \item{df2vec}{A named vector}
  \item{t2p2t}{List of species missing from phylo and vice versa}
  \item{internal2tips}{Vector of tips descended from a node}
  \item{node.age}{Phylo object with phylo\$ages vector of node ages corresponding to phylo\$edge}
  \item{phylo2phylog}{An ade4 phylog}
  \item{sortColumns}{A data.frame or matrix with columns sorted by name}
  \item{sortRows}{A data.frame or matrix with rows sorted by name}
  
}

\author{ Steven Kembel <skembel@berkeley.edu>, David Ackerly <dackerly@berkeley.edu> }
\keyword{manip}
