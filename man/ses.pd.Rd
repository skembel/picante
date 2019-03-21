\name{ses.pd}
\alias{ses.pd}

\title{ Standardized effect size of PD }
\description{
  Standardized effect size of phylogenetic diversity (Faith's PD) in communities.
}
\usage{
ses.pd(samp, tree, null.model = c("taxa.labels", "richness", "frequency",
    "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
    runs = 999, iterations = 1000, include.root=TRUE)
}

\arguments{
  \item{ samp }{ Community data matrix }
  \item{ tree }{ Phylogeny (phylo object) }
  \item{ null.model }{ Null model to use (see Details section for description) }
  \item{ runs }{ Number of randomizations }
  \item{ iterations }{ Number of iterations to use for each randomization (for independent swap and trial null models) }
  \item{ include.root }{ Include distance to root node in calculation of PD (see documentation in \code{\link{pd}} function }
}

\value{
  A data frame of results for each community
  \item{ntaxa}{Number of taxa in community}
  \item{pd.obs}{Observed PD in community}
  \item{pd.rand.mean}{Mean PD in null communities}
  \item{pd.rand.sd}{Standard deviation of PD in null communities}
  \item{pd.obs.rank}{Rank of observed PD vs. null communities}
  \item{pd.obs.z}{Standardized effect size of PD vs. null communities (= (pd.obs - pd.rand.mean) / pd.rand.sd)}
  \item{pd.obs.p}{P-value (quantile) of observed PD vs. null communities (= mpd.obs.rank / runs + 1)}
  \item{runs}{Number of randomizations}  
}
\details{
Currently implemented null models (arguments to null.model):
\describe{
    \item{taxa.labels}{ Shuffle taxa labels across tips of phylogeny (across all taxa included in phylogeny)}
    \item{richness}{ Randomize community data matrix abundances within samples (maintains sample species richness)}    
    \item{frequency}{ Randomize community data matrix abundances within species (maintains species occurence frequency)}    
    \item{sample.pool}{ Randomize community data matrix by drawing species from pool of species occurring in at least one community (sample pool) with equal probability}  
    \item{phylogeny.pool}{ Randomize community data matrix by drawing species from pool of species occurring in the phylogeny (phylogeny pool) with equal probability}  
    \item{independentswap}{ Randomize community data matrix with the independent swap algorithm (Gotelli 2000) maintaining species occurrence frequency and sample species richness }
    \item{trialswap}{ Randomize community data matrix with the trial-swap algorithm (Miklos & Podani 2004) maintaining species occurrence frequency and sample species richness }
    }
}

\references{ Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.

Proches, S., Wilson, J.R.U. and Cowling, R.M. 2006. How much evolutionary history in a 10 x 10m plot? Proceedings of Royal Society of London B, Biological Sciences 273:1143-1148.}
\author{ Steven Kembel <steve.kembel@gmail.com> }
\seealso{ \code{\link{pd}}, \code{\link{randomizeMatrix}}  }
\examples{
data(phylocom)
ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99)}
\keyword{univar}
