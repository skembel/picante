\name{randomizeSample}
\alias{randomizeSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Null models for community data matrix randomization}
\description{
  Various null models for randomizing community data matrices
}
\usage{
randomizeSample(samp, null.model = c("keepFreq", "keepRichness", "keepBoth"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{ Community data matrix }
  \item{null.model}{ Null model }
}
\value{
  Randomized community data matrix
}
\references{ Gotelli. 2000. Ecology. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\section{Warning }{keepBoth (independent swap) works ONLY with presence-absence}
\keyword{manip}
