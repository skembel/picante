\name{traitgram}
\alias{traitgram}

\title{ Draw phylogeny with nodes at trait positions }
\description{
  Draws a phylogeny where x position of nodes and tips corresponds to value of a continuous trait variable, and y position corresponds to node depth (i.e. age).
}
\usage{
traitgram(x, phy, xaxt = 's', underscore = FALSE, show.names = TRUE,
            show.xaxis.values = TRUE, method = c('ML','pic'), ...)
}

\arguments{
  \item{x}{ Trait vector (same order as phy\$tip.label, or with taxon names in names) }
  \item{phy}{ phylo object }
  \item{xaxt}{ x axis default type }
  \item{underscore}{ if FALSE remove underscore from taxonomic names }
  \item{show.names}{ if TRUE show taxon names across tips of phylogeny }
  \item{show.xaxis.values}{ if TRUE show values for trait on x=axis }
  \item{method}{ method for calculation of internal trait values. 'ML' = maximum likelihood method; 'pic' = independent contrasts method. pic option can be used when ML fails to converge or otherwise seems to fail to correctly reconstruct ancestral values}
  \item{...}{ Additional arguments passed to plot }
}
\value{
    Plots a traitgram, no values returned.
}
\author{ David Ackerly <dackerly@berkeley.edu> }
\references{ 
Ackerly, D. D. 2009. Conservatism and diversification of plant functional traits: Evolutionary rates versus phylogenetic signal. Proceedings of the National Academy of Sciences USA 106:19699-19706. doi: 10.1073/pnas.0901635106. 

Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue. 2009. Climate, Niche Evolution, and Diversification of the "bird-cage" Evening Primroses (Oenothera, Sections Anogra and Kleinia). American Naturalist 173:225-240. 
}
\examples{
randtree <- rcoal(20)
randtraits <- rTraitCont(randtree)
traitgram(randtraits,randtree)
traitgram(randtraits,randtree,method='pic')
}
\keyword{univar}