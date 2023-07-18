##' Phylogenetic index of beta-diversity PhyloSor
##' 
##' Fraction of branch-length shared between two communities
##' 
##' 
##' @param samp Community data matrix
##' @param tree Object of class phylo - a rooted phylogeny
##' @return A distance object of the PhyloSor index of similarity between
##' communities, the fraction of PD (branch-length) shared between two samples
##' @note The root of the supplied tree is included in calculations of
##' PhyloSor. The supplied tree must be rooted. Single-species samples will be
##' assigned a PD value equal to the distance from the root to the present.
##' @section Warning : The phylosor of all samples will include the branch
##' length connecting taxa in those samples and the root of the supplied tree.
##' The root of the supplied tree may not be spanned by any taxa in the sample.
##' If you want the root of your tree to correspond to the most recent ancestor
##' of the taxa actually present in your sample, you should prune the tree
##' before running \code{phylosor}:
##' 
##' \code{prunedTree <- prune.sample(sample,tree)}
##' @author Helene Morlon <morlon.helene@@gmail.com> and Steven Kembel
##' <steve.kembel@@gmail.com>
##' @seealso \code{\link{phylosor.rnd}}, \code{\link{pd}}
##' @references Bryant, J.B., Lamanna, C., Morlon, H., Kerkhoff, A.J., Enquist,
##' B.J., Green, J.L. 2008. Microbes on mountainsides: Contrasting elevational
##' patterns of bacterial and plant diversity. Proceedings of the National
##' Academy of Sciences 105 Supplement 1: 11505-11511
##' @keywords univar
##' @examples
##' 
##' data(phylocom)
##' phylosor(phylocom$sample, phylocom$phylo)
##' @export phylosor
phylosor <- function(samp, tree) {
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }

  if (!is.rooted(tree)) {
    stop("Rooted phylogeny required for phylosor calculation")
  }

  samp <- as.matrix(samp)
  s <- nrow(samp)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(samp)
  colnames(phylodist) <- rownames(samp)

  samp_comb <- matrix(NA, s * (s - 1) / 2, ncol(samp))
  colnames(samp_comb) <- colnames(samp)

  i <- 1
  for (l in 1:(s - 1))
  {
    for (k in (l + 1):s)
    {
      samp_comb[i, ] <- samp[l, ] + samp[k, ]
      i <- i + 1
    }
  }

  pdsamp <- pd(samp, tree)
  pdsamp_comb <- pd(samp_comb, tree)

  i <- 1
  for (l in 1:(s - 1)) {
    pdl <- pdsamp[l, "PD"]
    for (k in (l + 1):s) {
      pdk <- pdsamp[k, "PD"]
      pdcomb <- pdsamp_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
      phylodist[k, l] <- 2 * pdsharedlk / (pdl + pdk)
      i <- i + 1
    }
  }
  return(as.dist(phylodist))
}



##' Null PhyloSor values of phylogenetic beta-diversity
##' 
##' PhyloSor values obtained by randomization for different choices of null
##' models
##' 
##' Currently implemented null models (arguments to null.model): \describe{
##' \item{taxa.labels}{ Shuffle community data matrix labels. Maintains species
##' richness in each community and species shared between communities. Should
##' be used with cstSor=TRUE} \item{frequency}{ Randomize community data matrix
##' abundances within species (maintains species occurence frequency). Does not
##' maintain species richness in communities nor species shared between
##' communities. Can only be used with cstSor=FALSE} \item{richness}{ With
##' cstSor=TRUE: For each pair of community, maintains species richness in each
##' community and species shared between communities. Sample in the species
##' pool with equal probability; With cstSor=FALSE: Maintains species richness
##' in each community, does not maintain species shared between communities.
##' Sample in the species pool with equal probability} \item{independentswap}{
##' Randomize community data matrix with the independent swap algorithm
##' (Gotelli 2000) maintaining species occurrence frequency and sample species
##' richness. Can only be used with cstSor=FALSE} \item{trialswap}{ Randomize
##' community data matrix with the trial-swap algorithm (Miklos & Podani 2004)
##' maintaining species occurrence frequency and sample species richness. Can
##' only be used with cstSor=FALSE} }
##' 
##' @param samp Community data matrix
##' @param tree Object of class phylo - a rooted phylogeny
##' @param cstSor TRUE if the Sorensen similarity should be kept constant
##' across communities. FALSE otherwise
##' @param null.model Null model to use (see Details section)
##' @param runs Number of randomizations
##' @param iterations Number of iterations to use for each randomization (for
##' independent swap and trial null models)
##' @return A list of length the number of runs. Each element of the list is a
##' distance matrix containing the PhyloSor values of phylogenetic
##' beta-diversity obtained by randomization
##' @author Helene Morlon <morlon.helene@@gmail.com> and Steven Kembel
##' <steve.kembel@@gmail.com>
##' @seealso \code{\link{phylosor}}, \code{\link{randomizeMatrix}}
##' @references Bryant, J.B., Lamanna, C., Morlon, H., Kerkhoff, A.J., Enquist,
##' B.J., Green, J.L. 2008. Microbes on mountainsides: Contrasting elevational
##' patterns of bacterial and plant diversity. Proceedings of the National
##' Academy of Sciences 105 Supplement 1: 11505-11511
##' @keywords univar
##' @examples
##' 
##' data(phylocom)
##' phylosor.rnd(phylocom$sample,phylocom$phylo,cstSor=TRUE,null.model="richness",runs=5)
##' 
##' @export phylosor.rnd
phylosor.rnd <- function(samp, tree, cstSor = TRUE, null.model = c("taxa.labels", "frequency", "richness", "independentswap", "trialswap"), runs = 999, iterations = 1000) {
  Res <- list()

  if (cstSor == TRUE) {
    if (null.model == "taxa.labels") {
      for (r in 1:runs)
      {
        Res <- c(Res, list(.phylosor.taxaShuffle(samp, tree)))
      }
    } else if (null.model == "richness") {
      for (r in 1:runs)
      {
        Res <- c(Res, list(.phylosor.richness(samp, tree)))
      }
    } else {
      stop("This null model does not maintain Sorensen similarity: use cstSor=FALSE, or choose an other null model")
    }
  } else {
    if (null.model == "taxa.labels") {
      warning("This null model maintains Sorensen similarity")
      for (r in 1:runs)
      {
        Res <- c(Res, list(.phylosor.taxaShuffle(samp, tree)))
      }
    } else {
      for (r in 1:runs)
      {
        Res <- c(Res, list(phylosor(randomizeMatrix(samp, null.model), tree)))
      }
    }
  }

  return(Res)
}


##########################################################################################
.phylosor.taxaShuffle <- function(samp, tree) {
  sampr <- samp
  colnames(sampr) <- sample(colnames(samp))
  return(phylosor(sampr, tree))
}

##########################################################################################
.phylosor.richness <- function(samp, tree) {
  s <- nrow(samp)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(samp)
  colnames(phylodist) <- rownames(samp)

  for (l in 1:(s - 1))
  {
    for (k in (l + 1):s)
    {
      sampr <- samp
      colnames(sampr) <- sample(colnames(samp))
      pdl <- pd(sampr[l, , drop = FALSE], tree)$PD
      pdk <- pd(sampr[k, , drop = FALSE], tree)$PD
      pdtot <- pd((sampr[l, , drop = FALSE] + sampr[k, , drop = FALSE]), tree)$PD
      pdsharedlk <- pdl + pdk - pdtot
      phylodist[k, l] <- 2 * pdsharedlk / (pdl + pdk)
    }
  }
  return(as.dist(phylodist))
}
