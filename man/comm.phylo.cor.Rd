\name{comm.phylo.cor}
\alias{comm.phylo.cor}

\title{ Correlations between species co-occurrence and phylogenetic distances }
\description{
  Calculates measures of community phylogenetic structure (correlation between co-occurrence and phylogenetic distance) under various null models
}
\usage{
comm.phylo.cor(samp, phylo, metric = c("cij", "checkerboard","jaccard", "roij"),
    null.model = c("sample.taxa.labels", "pool.taxa.labels","frequency", "richness",
    "weighted.sample.pool"), runs = 99)

}
\arguments{
  \item{samp}{ Community data matrix }
  \item{phylo}{ Phylogenetic tree }
  \item{metric}{ Metric of co-occurrence to use }
  \item{null.model}{ Null model to use }
  \item{runs}{ Number of runs (randomizations) }
}
\value{
  \item{obs.corr }{Observed co-occurrence/phylogenetic distance correlation}
    obs.corr.p
    obs.rank
    runs = runs
    obs.rand.p
    random.corrs = vector(length = runs)    
  ...
}
\references{ Cavender-Bares et al. 2004, 2006 }
\author{ Steve Kembel <skembel@berkeley.edu> }
\section{Warning }{Weighted null currently only works with presence-absence}
\keyword{univar}
