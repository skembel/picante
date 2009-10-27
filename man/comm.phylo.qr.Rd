\name{comm.phylo.qr}
\alias{comm.phylo.qr}

\title{ Quantile regression slopes between species co-occurrence and phylogenetic distances }
\description{
  Calculates measures of community phylogenetic structure (quantile regression  between co-occurrence and phylogenetic distance) to patterns expected under various null models
}
\usage{
comm.phylo.qr(samp, phylo, metric = c("cij", "checkerboard", "jaccard", "doij"),
    null.model = c("sample.taxa.labels", "pool.taxa.labels",
        "frequency", "richness", "independentswap","trialswap"), 
        quant = 0.75, runs = 999, show.plot = FALSE, ...)
}
\arguments{
  \item{samp}{ Community data matrix }
  \item{phylo}{ Phylogenetic tree }
  \item{metric}{ Metric of co-occurrence to use (see \code{\link{species.dist}})}
  \item{null.model}{ Null model to use (see Details section for description) }
  \item{quant}{ Quantile of slope to be fit (using \code{\link[quantreg]{rq}}) }
  \item{runs}{ Number of runs (randomizations) }
  \item{show.plot}{ Option to display a plot of co-occurrence versus phylogenetic distance with quantile regression slope fit }
  \item{...}{ Additional arguments to \link{randomizeSample}}
}
\value{
    A list with elements:
    \item{obs.qr.intercept }{ Observed co-occurrence/phylogenetic distance quantile regression intercept}
    \item{obs.qr.slope }{ Observed co-occurrence/phylogenetic distance quantile regression slope}    
    \item{obs.qr.slope.p}{ P-value of observed quantile regression slope significance versus null model (calculated based on comparison with randomizations)}
    \item{obs.rank}{ Rank of observed quantile regression slope vs. random}
    \item{runs}{ Number of runs (randomizations) }
    \item{random.qr.slopes}{A vector of quantile regression slopes calculated for each randomization}
}
\details{
This function fits a quantile regression of co-occurrence versus phylogenetic distances separating species, and compares observed patterns to the patterns expected under some null model. The quantile regressions are fit using the \code{\link[quantreg]{rq}} function from the \code{quantreg} package.

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
\references{
Cavender-Bares J., D.A. Ackerly, D. Baum and F.A. Bazzaz. 2004. Phylogenetic overdispersion in Floridian oak communities, American Naturalist, 163(6):823-843.
Slingsby, J. A. and G. A. Verboom. 2006. Phylogenetic relatedness limits coexistence at fine spatial scales: evidence from the schoenoid sedges (Cyperaceae: Schoeneae) of the Cape Floristic Region, South Africa. The American Naturalist 168:14-27.
}
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{randomizeSample}} }
\examples{
data(phylocom)
comm.phylo.qr(phylocom$sample, phylocom$phylo, metric="cij",null.model="sample.taxa.labels")}
\keyword{univar}
