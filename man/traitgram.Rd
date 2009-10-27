\name{traitgram}
\alias{traitgram}

\title{ Draw phylogeny with nodes at trait positions }
\description{
  Draws a phylogeny where x position of nodes and tips corresponds to value of a continuous trait variable, and y position corresponds to node depth (i.e. age).
}
\usage{
traitgram(x, phy, xaxt = 's', underscore = FALSE, show.names = TRUE,
            show.xaxis.values = TRUE, method = c('ace','pic'), ...)
}

\arguments{
  \item{x}{ Trait vector (same order as phy\$tip.label, or with taxon names in names) }
  \item{phy}{ phylo object }
  \item{xaxt}{ x axis default type }
  \item{underscore}{ if FALSE remove underscore from taxonomic names }
  \item{show.names}{ if TRUE show taxon names across tips of phylogeny }
  \item{show.xaxis.values}{ if TRUE show values for trait on x=axis }
  \item{method}{ method for calculation of internal trait values. 'ace' = ace maximum likelihood method from ape; 'pic' = independent contrasts method. pic option can be used when ace fails to converge or otherwise seems to fail to correctly reconstruct ancestral values}
  \item{...}{ Additional arguments passed to plot }
}
\value{
    Plots a traitgram, no values returned.
}
\author{ David Ackerly <dackerly@berkeley.edu> }
\seealso{ \code{\link{pic3}} }
\examples{
randtree <- rcoal(20)
randtraits <- evolve.brownian(randtree)
traitgram(randtraits,randtree)
traitgram(randtraits,randtree,method='pic')
}
\keyword{univar}