\name{phylosor.rnd}
\alias{phylosor.rnd}

\title{ Null PhyloSor values of phylogenetic beta-diversity }
\description{
  PhyloSor values obtained by randomization for different choices of null models
}
\usage{
phylosor.rnd(samp,tree, cstSor=TRUE, null.model=c("taxa.labels",
    "frequency","richness","independentswap","trialswap"),
    runs=999, iterations=1000)
}

\arguments{
  \item{ samp }{ Community data matrix }
  \item{ tree }{ Object of class phylo}
  \item{ cstSor }{ TRUE if the Sorensen similarity should be kept constant across communities. FALSE otherwise }
  \item{ null.model }{ Null model to use (see Details section) }
  \item{ runs }{ Number of randomizations }
  \item{ iterations }{ Number of iterations to use for each randomization (for independent swap and trial null models) }
}

\value{
  A list of length the number of runs. Each element of the list is a distance matrix containing the PhyloSor values of phylogenetic beta-diversity obtained by randomization
}
\details{
Currently implemented null models (arguments to null.model):
\describe{
    \item{taxa.labels}{ Shuffle community data matrix labels. Maintains species richness in each community and species shared between communities. Should be used with cstSor=TRUE}
    \item{frequency}{ Randomize community data matrix abundances within species (maintains species occurence frequency). Does not maintain species richness in communities nor species shared between communities. Can only be used with cstSor=FALSE}
    \item{richness}{ With cstSor=TRUE: For each pair of community, maintains species richness in each community and species shared between communities. Sample in the species pool with equal probability;
    With cstSor=FALSE: Maintains species richness in each community, does not maintain species shared between communities. Sample in the species pool with equal probability}
    \item{independentswap}{ Randomize community data matrix with the independent swap algorithm (Gotelli 2000) maintaining species occurrence frequency and sample species richness. Can only be used with cstSor=FALSE}
    \item{trialswap}{ Randomize community data matrix with the trial-swap algorithm (Miklos & Podani 2004) maintaining species occurrence frequency and sample species richness. Can only be used with cstSor=FALSE}
    }
}
\references{Bryant, J.B., Lamanna, C., Morlon, H., Kerkhoff, A.J., Enquist, B.J., Green, J.L. 2008. Microbes on mountainsides: Contrasting elevational patterns of bacterial and plant diversity. Proceedings of the National Academy of Sciences 105 Supplement 1: 1505-11511}
\author{ Helene Morlon <morlon.helene@gmail.com> and Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{phylosor}}, \code{\link{randomizeSample}}  }
\examples{
data(phylocom)
phylosor.rnd(phylocom$sample,phylocom$phylo,cstSor=TRUE,null.model="richness",runs=5)
}
\keyword{univar}
