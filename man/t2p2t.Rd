\name{t2p2t}
\alias{t2p2t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ List of non-matching species from phylo and trait data objects }
\description{
  List of non-matching species from phylo and trait data objects
}
\usage{
t2p2t(phy, species)
}

\arguments{
  \item{phy}{ phylo object }
  \item{species}{ Vector of species names }
}

\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{missing.from.phy:}{Species in trait but not phylo}
  \item{missing.from.species: }{Species in phylo but not trait}
  ...
}

\author{ David Ackerly <dackerly@berkeley.edu> }
\keyword{manip}
