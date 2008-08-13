\name{randomizeSample}
\alias{randomizeSample}

\title{ Null models for community data matrix randomization}
\description{
  Various null models for randomizing community data matrices
}
\usage{
randomizeSample(samp, null.model = c("frequency", "richness", "independentswap","trialswap"), lang=c("C","R"), it=1000)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{null.model}{ Null model
    \item{frequency}{Randomize community data matrix abundances within species (maintains species occurence frequency)}
    \item{richness}{Randomize community data matrix abundances within samples (maintains sample species richness)}
    \item{independentswap}{Randomize community data matrix with the independent-swap algorithm (Gotelli 2000) by drawing species from sample pool with probability weighted by occurrence frequency (maintains both species occurrence frequency and sample species richness)}
    \item{trialswap}{Randomize community data matrix with the trial-swap algorithm (Miklos & Podani 2004) by drawing species from sample pool with probability weighted by occurrence frequency (maintains both species occurrence frequency and sample species richness)}
  }
  \item{lang}{run code in C or R? C is faster.} 
  \item{it}{ number of independent or trial-swaps to perform }
}
\value{
  Randomized community data matrix
}
\references{Gotelli, N.J. 2000. Null model analysis of species co-occurrence patterns. Ecology 81: 2606-2621.
\cr
\cr
Miklos I. & Podani J. 2004. Randomization of presence-absence matrices: Comments and new algorithms. Ecology 85: 86-92.
}
\author{ Steve Kembel <skembel@berkeley.edu> }
\examples{
data(phylocom)
randomizeSample(phylocom$sample, null.model="richness")}
\keyword{manip}
