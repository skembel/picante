##' Rao's quadratic entropy
##' 
##' Calculates Rao's quadratic entropy, a measure of within- and
##' among-community diversity taking species dissimilarities into account
##' 
##' Rao's quadratic entropy (Rao 1982) is a measure of diversity in ecological
##' communities that can optionally take species differences (e.g. phylogenetic
##' dissimilarity) into account. This method is conceptually similar to
##' analyses of genetic diversity among populations (Nei 1973), but instead of
##' diversity of alleles among populations, it measures diversity of species
##' among communities.
##' 
##' If no phylogeny is supplied, Dkk is equivalent to Simpson's diversity
##' (probability that two individuals drawn from a community are from different
##' taxa), Dkl is a beta-diversity equivalent of Simpson's diversity
##' (probability that individuals drawn from each of two communities belong to
##' different taxa), and H is Dkl standardized to account for within-community
##' diversity.
##' 
##' If an ultrametric phylogeny is supplied, Dkk is equivalent to the mean
##' pairwise phylogenetic distance (distance to MRCA) between two individuals
##' drawn from a community, Dkl is the mean pairwise phylogenetic distance
##' between individuals drawn from each of two communities, and H is Dkl
##' standardized to account for within-community diversity.
##' 
##' \deqn{Dkl = sum(tij * xki * xli)}
##'
##' where \emph{xki} is the relative abundance of taxon \emph{i} in community
##' \emph{k} and \emph{tij} is a matrix of weights for all pairs of taxa
##' \emph{i,j}. Without a phylogeny, when \emph{i=j}, \emph{tij=0}, otherwise
##' \emph{tij=1}. With a phylogeny, \emph{tij} is the phylogenetic distance
##' to MRCA for taxa \emph{i,j}.
##' 
##' 
##' \deqn{Hkl = Dkl-(Dkk + Dll)/2}
##' 
##' Alpha, beta and total measure the average diversity within, among, and
##' across all communities based on Dkk and H values taking variation in number
##' of individuals per community into account. A Fst-like measure is calculated
##' by dividing beta by the total diversity across all samples.
##' 
##' @param comm Community data matrix
##' @param phy Object of class phylo - an ultrametric phylogenetic tree
##' (optional)
##' @return A list of results \item{ Dkk }{ Within-community diversity } \item{
##' Dkl }{ Among-community diversity } \item{ H }{ Among-community diversities
##' excluding within-community diversity } \item{ total }{ Total diversity
##' across all samples } \item{ alpha }{ Alpha diversity - average
##' within-community diversity } \item{ beta }{ Beta diversity - average
##' among-community diversity } \item{ Fst }{ Beta diversity / total diversity
##' }
##' @section Warning : Alpha, beta, and total diversity components and Fst
##' should not be interpreted as a measure of relative differentiation among
##' versus within communities. See Jost (2007) for a detailed description of
##' this problem. Hardy and Jost (2008) suggest Fst can be interpreted as
##' 'local species identity excess' or 'local phylogenetic similarity excess'
##' rather than as a measure of among-community differentiation.
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @seealso \code{\link{mpd}}, \code{\link{comdist}}
##' @references Hardy, O.J., and Jost. L. 2008. Interpreting and estimating
##' measures of community phylogenetic structuring. J. Ecol. 96:849-852.
##' 
##' Jost, L. 2007. Partitioning diversity into independent alpha and beta
##' components. Ecology 88: 24272439.
##' 
##' Nei, M. 1973. Analysis of gene diversity in sub-divided populations.
##' Proceedings of the National Academy of Sciences of the USA 70:3321-3323.
##' 
##' Rao, C.R. 1982. Diversity and dissimilarity coefficients: a unified
##' approach. Theoretical Population Biology 21:2443.
##' 
##' Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for
##' the analysis of phylogenetic community structure and trait evolution.
##' Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
##' @keywords univar
##' @examples
##' 
##' data(phylocom)
##' raoD(phylocom$sample)
##' raoD(phylocom$sample, phylocom$phylo)
##' 
##' @export raoD
raoD <- function(comm, phy = NULL) {
  res <- list()

  if (is.null(phy)) {
    tij <- 1 - diag(x = rep(1, length(comm[1, ])))
  } else {
    if (!is.ultrametric(phy)) {
      stop("Phylogeny must be ultrametric")
    }
    dat <- match.phylo.comm(phy, comm)
    comm <- dat$comm
    phy <- dat$phy
    tij <- cophenetic(phy) / 2
  }

  x <- as.matrix(comm)
  S <- length(x[1, ])
  N <- length(x[, 1])
  total <- apply(x, 1, sum)
  samp.relabund <- total / sum(x)
  x.combined <- matrix(apply(x, 2, sum), nrow = 1) / sum(x)
  x <- sweep(x, 1, total, "/")

  D <- vector(length = N)
  names(D) <- rownames(x)
  for (k in 1:N) {
    D[k] <-
      sum(tij * outer(as.vector(t(x[k, ])), as.vector(t(x[k, ]))))
  }
  res$Dkk <- D

  Dkl <- matrix(nrow = N, ncol = N)
  for (k in 1:N) {
    for (l in 1:N) {
      Dkl[k, l] <-
        sum(tij * outer(as.vector(t(x[k, ])), as.vector(t(x[l, ]))))
    }
  }

  row.names(Dkl) <- row.names(x)
  colnames(Dkl) <- row.names(x)
  H <- Dkl
  res$Dkl <- Dkl

  for (k in 1:N) {
    for (l in 1:N) {
      H[k, l] <- Dkl[k, l] - (Dkl[k, k] + Dkl[l, l]) / 2
    }
  }
  res$H <- H

  res$total <- sum(tij * outer(as.vector(t(x.combined)), as.vector(t(x.combined))))

  res$alpha <- sum(res$Dkk * samp.relabund)

  res$beta <- res$total - res$alpha

  res$Fst <- res$beta / res$total

  return(res)
}
