\name{comm.phylo.cor}
\alias{comm.phylo.cor}

\title{ Correlations between species co-occurrence and phylogenetic distances }
\description{
  Calculates measures of community phylogenetic structure (correlation between co-occurrence and phylogenetic distance) to patterns expected under various null models
}
\usage{
comm.phylo.cor(samp, phylo, metric = c("cij", "checkerboard", "jaccard", "doij"),
    null.model = c("sample.taxa.labels", "pool.taxa.labels",
        "frequency", "richness", "independentswap","trialswap"), runs = 999, ...)
}
\arguments{
  \item{samp}{ Community data matrix }
  \item{phylo}{ Phylogenetic tree }
  \item{metric}{ Metric of co-occurrence to use (see \code{\link{species.dist}})}
  \item{null.model}{ Null model to use (see Details section for description) }
  \item{runs}{ Number of runs (randomizations) }
  \item{...}{ Additional arguments to \link{randomizeSample}}
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
\details{
Currently implemented null models (arguments to null.model):
\describe{
        \item{sample.taxa.labels}{Shuffle phylogeny tip labels (only within set of taxa present in community data)}
        \item{pool.taxa.labels}{Shuffle phylogeny tip labels (across all taxa included in phylogenetic tree)}
        \item{frequency}{Randomize community data matrix abundances within species (maintains species occurence frequency)}
        \item{richness}{Randomize community data matrix abundances within samples (maintains sample species richness)}
        \item{independentswap}{Randomize community data matrix maintaining species occurrence frequency and site richnessing using independent swap}
        \item{trialswap}{Randomize community data matrix maintaining species occurrence frequency and site richnessing using trial swap}     
    }
}
\references{ Cavender-Bares J., D.A. Ackerly, D. Baum and F.A. Bazzaz. 2004. Phylogenetic overdispersion in Floridian oak communities, American Naturalist, 163(6):823-843. }
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{randomizeSample}} }
\examples{
data(phylocom)
comm.phylo.cor(phylocom$sample, phylocom$phylo, metric="cij",null.model="sample.taxa.labels")}
\keyword{univar}
