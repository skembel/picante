\name{raoD}
\alias{raoD}

\title{ Rao's quadratic entropy }
\description{
  Calculates Rao's quadratic entropy, a measure of within- and among-community diversity taking species dissimilarities into account
}
\usage{
  raoD(comm, phy=NULL)
}

\arguments{
  \item{ comm }{ Community data matrix }
  \item{ phy }{ Object of class phylo - an ultrametric phylogenetic tree (optional) }
}

\value{
  A list of results
  \item{ Dkk }{ Within-community diversity }
  \item{ Dkl }{ Among-community diversity }
  \item{ H }{ Among-community diversities excluding within-community diversity }
  \item{ total }{ Total diversity across all samples }
  \item{ alpha }{ Alpha diversity - average within-community diversity }
  \item{ beta }{ Beta diversity - average among-community diversity }
  \item{ Fst }{ Beta diversity / total diversity }
}

\details{
Rao's quadratic entropy (Rao 1982) is a measure of diversity in ecological communities that can optionally take species differences (e.g. phylogenetic dissimilarity) into account. This method is conceptually similar to analyses of genetic diversity among populations (Nei 1973), but instead of diversity of alleles among populations, it measures diversity of species among communities. 

If no phylogeny is supplied, Dkk is equivalent to Simpson's diversity (probability that two individuals drawn from a community are from different taxa), Dkl is a beta-diversity equivalent of Simpson's diversity (probability that individuals drawn from each of two communities belong to different taxa), and H is Dkl standardized to account for within-community diversity.

If an ultrametric phylogeny is supplied, Dkk is equivalent to the mean pairwise phylogenetic distance (distance to MRCA) between two individuals drawn from a community, Dkl is the mean pairwise phylogenetic distance between individuals drawn from each of two communities, and H is Dkl standardized to account for within-community diversity.

\emph{D[kl] = sum(t[ij] * x[ki] * x[lj])}
    
where \emph{x[ki]} is the relative abundance of taxon \emph{i} in community \emph{k} and \emph{t[ij]} is a matrix of weights for all pairs of taxa \emph{i,j}. Without a phylogeny, when \emph{i=j}, \emph{t[ij]=0}, otherwise \emph{t[ij]=1}. With a phylogeny, \emph{t[ij]} is the phylogenetic distance to MRCA for taxa \emph{i,j}.
            
\emph{H[kl] = D[kl] - (D[kk] + D[ll])/2}

Alpha, beta and total measure the average diversity within, among, and across all communities based on Dkk and H values taking variation in number of individuals per community into account. A Fst-like measure is calculated by dividing beta by the total diversity across all samples.
}

\section{ Warning }{
Alpha, beta, and total diversity components and Fst should not be interpreted as a measure of relative differentiation among versus within communities. See Jost (2007) for a detailed description of this problem. Hardy and Jost (2008) suggest Fst can be interpreted as 'local species identity excess' or 'local phylogenetic similarity excess' rather than as a measure of among-community differentiation.
}
\references{
Hardy, O.J., and Jost. L. 2008. Interpreting and estimating measures of community phylogenetic structuring. J. Ecol. 96:849-852.

Jost, L. 2007. Partitioning diversity into independent alpha and beta components. Ecology 88: 24272439.

Nei, M. 1973. Analysis of gene diversity in sub-divided populations. Proceedings of the National Academy of Sciences of the USA 70:3321-3323.

Rao, C.R. 1982. Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology 21:2443.

Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
}
\seealso{ \code{\link{mpd}}, \code{\link{comdist}}  }
\author{ Steven Kembel <skembel@uoregon.edu> }
\examples{
data(phylocom)
raoD(phylocom$sample)
raoD(phylocom$sample, phylocom$phylo)
}
\keyword{univar}