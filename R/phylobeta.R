# phylobeta diversity functions

# comdist
# mean pairwise distance among taxa from two communities


##' Calculates inter-community mean pairwise distance
##' 
##' Calculates MPD (mean pairwise distance) separating taxa in two communities,
##' a measure of phylogenetic beta diversity
##' 
##' This function calculates a measure of phylogenetic beta diversity: the
##' expected phylogenetic distance separating two individuals or taxa drawn
##' randomly from different communities.
##' 
##' @param comm Community data matrix
##' @param dis Interspecific distance matrix
##' @param abundance.weighted Should mean pairwise distances separating species
##' in two communities be weighted by species abundances? (default = FALSE)
##' @param threads is an interger of number of threads for parallel processing
##' @return Distance object of MPD values separating each pair of communities.
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @seealso \code{\link{mpd}}, \code{\link{ses.mpd}}
##' @references C.O. Webb, D.D. Ackerly, and S.W. Kembel. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Bioinformatics 18:2098-2100.
##' @keywords univar
##' @examples
##' 
##' data(phylocom)
##' comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
##' @export comdist
comdist <- function(comm, dis, abundance.weighted = FALSE, Rcpp = TRUE, threads = 1) {
  dat <- match.comm.dist(comm, dis)
  x <- dat$comm
  dis <- as.matrix(dat$dist)

  if (!abundance.weighted) {
    x <- decostand(x, method = "pa")
  }
  N <- dim(x)[1]
  x <- decostand(x, method = "total", MARGIN = 1)


if(Rcpp){
  comdist <- matrix(nrow = N, ncol = N)
  for (l in 1:(N - 1)) {
    for (k in (l + 1):N) {
      comdist[k, l] <-
        sum(dis * outer(as.vector(t(x[k, ])), as.vector(t(x[l, ]))))
    }
  }

} else {
  comdist <- distmat_rcpp(as.matrix(x), dis, threads = threads)
}

  row.names(comdist) <- row.names(x)
  colnames(comdist) <- row.names(x)
  return(as.dist(comdist))
}


# comdistnt
# mean distance to closest relative between taxa from two communities


##' Calculates inter-community mean nearest taxon distance
##' 
##' Calculates MNTD (mean nearest taxon distance) separating taxa in two
##' communities, a measure of phylogenetic beta diversity
##' 
##' This metric has also been referred to as MNND (mean nearest neighbour
##' distance).
##' 
##' This function calculates a measure of phylogenetic beta diversity: the
##' average phylogenetic distance to the most similar taxon or individual in
##' the other community for taxa or individuals in two communities.
##' 
##' @aliases comdistnt comdistnn
##' @param comm Community data matrix
##' @param dis Interspecific distance matrix
##' @param abundance.weighted Should mean nearest taxon distances from each
##' species to species in the other community be weighted by species abundance?
##' (default = FALSE)
##' @param exclude.conspecifics Should conspecific taxa in different
##' communities be exclude from MNTD calculations? (default = FALSE)
##' @return Distance object of MNTD values separating each pair of communities.
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @seealso \code{\link{mntd}}, \code{\link{ses.mntd}}
##' @references C.O. Webb, D.D. Ackerly, and S.W. Kembel. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Bioinformatics 18:2098-2100.
##' @keywords univar
##' @examples
##' 
##' data(phylocom)
##' comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)
##' @export comdistnt
comdistnt <- function(comm, dis, abundance.weighted = FALSE, exclude.conspecifics = FALSE) {
  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]
  comm <- decostand(comm, method = "total", MARGIN = 1)
  comdisnt <- matrix(nrow = N, ncol = N)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      sppInSample1 <- colnames(comm[i, comm[i, ] > 0, drop = FALSE])
      sppInSample2 <- colnames(comm[j, comm[j, ] > 0, drop = FALSE])
      if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 1)) {
        sample.dis <- dis[sppInSample1, sppInSample2, drop = FALSE]
        if (exclude.conspecifics) {
          sample.dis[sample.dis == 0] <- NA
        }
        # TODO fix min throws errors on empty set
        sample1NT <- apply(sample.dis, 1, min, na.rm = TRUE)
        sample1NT[sample1NT == Inf] <- NA
        sample2NT <- apply(sample.dis, 2, min, na.rm = TRUE)
        sample2NT[sample2NT == Inf] <- NA
        if (abundance.weighted) {
          sample1.weights <- as.numeric(comm[i, sppInSample1])
          sample2.weights <- as.numeric(comm[j, sppInSample2])
          if (any(is.na(sample1NT))) {
            miss <- which(is.na(sample1NT))
            sample1NT <- sample1NT[-miss]
            sample1.weights <- sample1.weights[-miss]
            sample1.weights <- sample1.weights / sum(sample1.weights)
          }
          if (any(is.na(sample2NT))) {
            miss <- which(is.na(sample2NT))
            sample2NT <- sample2NT[-miss]
            sample2.weights <- sample2.weights[-miss]
            sample2.weights <- sample2.weights / sum(sample2.weights)
          }
          sampleNT <- c(sample1NT, sample2NT)
          sample.weights <- c(sample1.weights, sample2.weights)
          comdisnt[i, j] <- weighted.mean(sampleNT, sample.weights, na.rm = TRUE)
        } else {
          comdisnt[i, j] <- mean(c(sample1NT, sample2NT), na.rm = TRUE)
        }
      } else {
        comdisnt[i, j] <- NA
      }
    }
  }
  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  return(as.dist(t(comdisnt)))
}

# unifrac


##' Unweighted UniFrac distance between communities
##' 
##' Calculates unweighted UniFrac, a phylogenetic beta diversity metric of the
##' the unique (non-shared) fraction of total phylogenetic diversity
##' (branch-length) between two communities.
##' 
##' 
##' @param comm Community data matrix
##' @param tree Object of class phylo - a rooted phylogeny
##' @return A dist object of the unweighted UniFrac distances between
##' communities (the unique (non-shared) fraction of total phylogenetic
##' diversity (branch-length) between two communities).
##' @note The supplied tree must be rooted. Single-species samples will be
##' assigned a PD value equal to the distance from the root to the present.
##' @section Warning : The UniFrac distance between samples will include the
##' branch length connecting taxa in those samples and the root of the supplied
##' tree. The root of the supplied tree may not be spanned by any taxa in the
##' sample. If you want the root of your tree to correspond to the most recent
##' ancestor of the taxa actually present in your samples, you should prune the
##' tree before running \code{unifrac}: \code{prunedTree <-
##' prune.sample(sample,tree)}
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @seealso \code{\link{pd}}
##' @references Lozupone, C., Hamady, M., and Knight, R. 2006. UniFrac - an
##' online tool for comparing microbial community diversity in a phylogenetic
##' context. BMC Bioinformatics 7:371.
##' @keywords univar
##' @examples
##' 
##' data(phylocom)
##' unifrac(phylocom$sample, phylocom$phylo)
##' @export unifrac
unifrac <- function(comm, tree) {
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute UniFrac")
  }

  if (!is.rooted(tree)) {
    stop("Rooted phylogeny required for UniFrac calculation")
  }

  comm <- as.matrix(comm)
  s <- nrow(comm)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(comm)
  colnames(phylodist) <- rownames(comm)

  comm_comb <- matrix(NA, s * (s - 1) / 2, ncol(comm))
  colnames(comm_comb) <- colnames(comm)

  i <- 1
  for (l in 1:(s - 1))
  {
    for (k in (l + 1):s)
    {
      comm_comb[i, ] <- comm[l, ] + comm[k, ]
      i <- i + 1
    }
  }

  pdcomm <- pd(comm, tree)
  pdcomm_comb <- pd(comm_comb, tree)

  i <- 1
  for (l in 1:(s - 1)) {
    pdl <- pdcomm[l, "PD"]
    for (k in (l + 1):s) {
      pdk <- pdcomm[k, "PD"]
      pdcomb <- pdcomm_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
      phylodist[k, l] <- (pdcomb - pdsharedlk) / pdcomb
      i <- i + 1
    }
  }
  return(as.dist(phylodist))
}
