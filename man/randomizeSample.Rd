\name{randomizeSample}
\alias{randomizeSample}

\title{ Null models for community data matrix randomization}
\description{
  Various null models for randomizing community data matrices
}
\usage{
randomizeSample(samp, null.model = c("frequency", "richness", "both"))
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{null.model}{ Null model
    \item{frequency}{Randomize community data matrix abundances within species (maintains species occurence frequency)}
    \item{richness}{Randomize community data matrix abundances within samples (maintains sample species richness)}
    \item{both}{Randomize community data matrix by drawing species from sample pool with probability weighted by occurrence frequency (maintains both species occurrence frequency and sample species richness)}
  }
}
\value{
  Randomized community data matrix
}
\references{Gotelli, N.J. 2000. Null model analysis of species co-occurrence patterns. Ecology 81: 2606-2621.}
\author{ Steve Kembel <skembel@berkeley.edu> }
\section{Warning }{ Null model \code{both} currently only works with presence-absence data. Convert your data to presence-absence before using this null model (e.g. \code{decostand(x,method="pa")})}
\examples{
data(phylocom)
randomizeSample(phylocom$sample, null.model="richness")}
\keyword{manip}
