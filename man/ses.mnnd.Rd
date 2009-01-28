\name{ses.mnnd}
\alias{ses.mnnd}

\title{ Standardized effect size of MNND }
\description{
  Standardized effect size of mean nearest neighbour distances in communities. When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
}
\usage{
ses.mnnd(samp, dis, null.model = c("taxa.labels","sample.pool", "phylogeny.pool",
            "independentswap", "trialswap"), runs = 999, iterations = 1000)    
}

\arguments{
  \item{ samp }{ Community data matrix }
  \item{ dis }{ Distance matrix (generally a phylogenetic distance matrix)}
  \item{ null.model }{ Null model to use (see Details section for description) }
  \item{ runs }{ Number of randomizations }
  \item{ iterations }{ Number of iterations to use for each randomization (for independent swap and trial null models) }
}

\value{
  A data frame of results for each community
  \item{ntaxa}{Number of taxa in community}
  \item{mnnd.obs}{Observed MNND in community}
  \item{mnnd.rand.mean}{Mean MNND in null communities}
  \item{mnnd.rand.sd}{Standard deviation of MNND in null communities}
  \item{mnnd.obs.rank}{Rank of observed MNND vs. null communities}
  \item{mnnd.obs.z}{Standardized effect size of MNND vs. null communities (= (mnnd.obs - mnnd.rand.mean) / mnnd.rand.sd, equivalent to -NTI)}
  \item{mnnd.obs.p}{P-value (quantile) of observed MNND vs. null communities (= mnnd.obs.rank / runs + 1)}
  \item{runs}{Number of randomizations}  
}
\details{
Currently implemented null models (arguments to null.model):
\describe{
    \item{taxa.labels}{ Shuffle distance matrix labels (across all taxa included in distance matrix)}
    \item{sample.pool}{ Randomize community data matrix by drawing species from pool of species occurring in at least one community (sample pool) with equal probability}  
    \item{phylogeny.pool}{ Randomize community data matrix by drawing species from pool of species occurring in the distance matrix (phylogeny pool) with equal probability}  
    \item{frequency}{ Randomize community data matrix abundances within species (maintains species occurence frequency)}
    \item{richness}{ Randomize community data matrix abundances within samples (maintains sample species richness)}
    \item{independentswap}{ Randomize community data matrix with the independent swap algorithm (Gotelli 2000) maintaining species occurrence frequency and sample species richness }
    \item{trialswap}{ Randomize community data matrix with the trial-swap algorithm (Miklos & Podani 2004) maintaining species occurrence frequency and sample species richness }
}
}
\references{ Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}. }
\author{ Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{mnnd}}, \code{\link{randomizeSample}}  }
\examples{
data(phylocom)
ses.mnnd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels")}
\keyword{univar}
