\name{taxaShuffle}
\alias{taxaShuffle}

\title{ Shuffle taxa names on a dist or matrix object}
\description{
  Shuffle taxa names on a dist or matrix object
}
\usage{
taxaShuffle(x, strata)
}

\arguments{
  \item{x}{ Matrix or dist object of interspecific distances }
  \item{strata}{ Strata for randomization }
}
\details{
  Code from \code{\link[vegan]{mantel}}
}
\value{
  Matrix with taxa names shuffled
}

\seealso{ \code{\link[vegan]{mantel}} }

\author{ Steve Kembel <skembel@berkeley.edu> }

\keyword{manip}
