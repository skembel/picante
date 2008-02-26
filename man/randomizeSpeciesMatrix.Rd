\name{randomizeSpeciesMatrix}
\alias{randomizeSpeciesMatrix}

\title{ Randomize community matrix }
\description{
  Randomize community matrix keeping sample richness constant and with option to keep species frequency constant
}
\usage{
randomizeSpeciesMatrix(x, keepSppFreq = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Community data matrix }
  \item{keepSppFreq}{ Should species frequencies be maintained? (TRUE or FALSE) }
}

\value{
  Randomized community data matrix
}
\references{ Gotelli. 2000. Ecology. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\section{Warning }{Slow! Works only with presence-absence matrices.}
\keyword{manip}


