\name{comm.phylo.cor}
\alias{comm.phylo.cor}

\title{ Correlations between species co-occurrence and phylogenetic distances }
\description{
  Calculates measures of community phylogenetic structure (correlation between co-occurrence and phylogenetic distance) to patterns expected under various null models
}
\usage{
comm.phylo.cor(samp, phylo,
    metric = c("cij", "checkerboard","jaccard", "roij"),
    null.model = c("sample.taxa.labels", "pool.taxa.labels","frequency",
    "richness", "weighted.sample.pool"), runs = 99)

}
\arguments{
  \item{samp}{ Community data matrix }
  \item{phylo}{ Phylogenetic tree }
  \item{metric}{ Metric of co-occurrence to use (see \code{\link{species.dist}})}
  \item{null.model}{ Null model to use:
        \item{sample.taxa.labels}{Shuffle phylogeny tip labels (only within set of taxa present in community data)}
        \item{pool.taxa.labels}{Shuffle phylogeny tip labels (across all taxa included in phylogenetic tree)}
        \item{frequency}{Randomize community data matrix abundances within species (maintains species occurence frequency)}
        \item{richness}{Randomize community data matrix abundances within samples (maintains sample species richness)}
        \item{weighted.sample.pool}{Randomize community data matrix by drawing species from sample pool with probability weighted by occurrence frequency}
  }
  \item{runs}{ Number of runs (randomizations) }
}
\value{
    A list with elements:
    \item{obs.corr }{ Observed co-occurrence/phylogenetic distance correlation}
    \item{obs.corr.p}{ P-value of observed correlation (standard P-value for correlation coefficient, not based on comparison with randomizations)}
    \item{obs.rank}{ Rank of observed correlation vs. random}
    \item{runs}{ Number of runs (randomizations) }
    \item{obs.rand.p}{ P-value of observed correlation vs. randomizations (= obs.rank / (runs + 1))}
    \item{random.corrs}{A vector of random correlation calculated for each run}
}
\references{ Cavender-Bares J., D.A. Ackerly, D. Baum and F.A. Bazzaz. 2004. Phylogenetic overdispersion in Floridian oak communities, American Naturalist, 163(6):823-843. }
\author{ Steve Kembel <skembel@berkeley.edu> }
\section{Warning }{ Weighted null model currently only works with presence-absence data. Convert your data to presence-absence before using this null model (e.g. \code{decostand(x,method="pa")})}
\examples{
data(phylocom)
comm.phylo.cor(phylocom$sample, phylocom$phylo, metric="cij",null.model="sample.taxa.labels")}
\keyword{univar}
