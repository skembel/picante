\name{ses.mntd}
\alias{ses.mntd}
\alias{ses.mnnd}

\title{ Standardized effect size of MNTD }
\description{
  Standardized effect size of mean nearest taxon distances in communities. When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
}
\usage{
ses.mntd(samp, dis, null.model = c("taxa.labels", "richness", "frequency",
    "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
    abundance.weighted=FALSE, runs = 999, iterations = 1000)
}

\arguments{
  \item{ samp }{ Community data matrix }
  \item{ dis }{ Distance matrix (generally a phylogenetic distance matrix)}
  \item{ null.model }{ Null model to use (see Details section for description) }
  \item{ abundance.weighted }{ Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE) }
  \item{ runs }{ Number of randomizations }
  \item{ iterations }{ Number of iterations to use for each randomization (for independent swap and trial null models) }
}

\value{
  A data frame of results for each community
  \item{ntaxa}{Number of taxa in community}
  \item{mntd.obs}{Observed MNTD in community}
  \item{mntd.rand.mean}{Mean MNTD in null communities}
  \item{mntd.rand.sd}{Standard deviation of MNTD in null communities}
  \item{mntd.obs.rank}{Rank of observed MNTD vs. null communities}
  \item{mntd.obs.z}{Standardized effect size of MNTD vs. null communities (= (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NTI)}
  \item{mntd.obs.p}{P-value (quantile) of observed MNTD vs. null communities (= mntd.obs.rank / runs + 1)}
  \item{runs}{Number of randomizations}  
}
\details{
The metric used by this function has also been referred to as MNND (mean nearest neighbour distance), and the function was named \code{ses.mnnd} in picante versions < 0.7.

Currently implemented null models (arguments to null.model):
\describe{
    \item{taxa.labels}{ Shuffle distance matrix labels (across all taxa included in distance matrix)}
    \item{richness}{ Randomize community data matrix abundances within samples (maintains sample species richness)}    
    \item{frequency}{ Randomize community data matrix abundances within species (maintains species occurence frequency)}
    \item{sample.pool}{ Randomize community data matrix by drawing species from pool of species occurring in at least one community (sample pool) with equal probability}  
    \item{phylogeny.pool}{ Randomize community data matrix by drawing species from pool of species occurring in the distance matrix (phylogeny pool) with equal probability}  
    \item{independentswap}{ Randomize community data matrix with the independent swap algorithm (Gotelli 2000) maintaining species occurrence frequency and sample species richness }
    \item{trialswap}{ Randomize community data matrix with the trial-swap algorithm (Miklos & Podani 2004) maintaining species occurrence frequency and sample species richness }
}
}
\references{ Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}. }
\author{ Steven Kembel <steve.kembel@gmail.com> }
\seealso{ \code{\link{mntd}}, \code{\link{randomizeMatrix}}  }
\examples{
data(phylocom)
ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels")}
\keyword{univar}
