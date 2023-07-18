# phyEstimate.R
# use phylogeny to predict trait values for new species


##' Phylogenetic estimation of traits for unobserved taxa
##'
##' Uses phylogenetic ancestral state reconstruction to estimate trait values
##' for unobserved taxa.
##'
##' These functions use phylogenetic ancestral state estimation to infer trait
##' values for novel taxa on a phylogenetic tree, for continuous
##' (\code{phyEstimate}) and discrete (\code{phyEstimateDisc}) traits.
##'
##' The required input is a phylogenetic tree object plus a vector or
##' data.frame containing estimated trait values for a subset of the taxa in
##' the phylogenetic tree. Trait values for taxa that are present in the tree
##' but not the trait data will be estimated using ancestral state estimation
##' (Garland and Ives 2000). Briefly, for each taxon present in the tree but
##' not the trait data, the phylogeny is rerooted at the most recent common
##' ancestor of the novel taxon and the rest of the phylogeny, and the trait
##' value of the novel taxon is estimated from the reconstructed trait value at
##' the root of the rerooted phylogeny.
##'
##' For \code{phyEstimateDisc}, the state with the highest support will be
##' reported if argument \code{best.state=TRUE}. If the best-supported state's
##' support is less than the specified \code{cutoff}, no best state is reported
##' and a \code{NA} value will be returned.
##'
##' @aliases phyEstimate phyEstimateDisc
##' @param phy phylo object
##' @param trait vector or data.frame containing trait values
##' @param method ancestral state estimation method used by \code{ace}
##' (default="pic")
##' @param best.state estimate best-supported trait state for discrete
##' variables? (default=TRUE)
##' @param cutoff support cutoff required to declare a best.state
##' @param ...  Additional arguments passed to \code{ace}
##' @return phyEstimate produces a data frame with columns: \item{est}{
##' Estimated trait value } \item{se}{ Standard error of estimated trait value
##' } phyEstimateDisc produces a data frame with columns: \item{states 1..N}{ A
##' column with statistical support is produced for each discrete trait state }
##' \item{estimated.state}{ If best.state=TRUE, a column with the state with
##' the highest support } \item{estimated.state.support}{ Statistical support
##' for the state with the highest support }
##' @author Steven Kembel <steve.kembel@gmail.com>
##' @references T. Garland Jr., and A.R. Ives. 2000. Using the past to predict
##' the present: confidence intervals for regression equations in phylogenetic
##' comparative methods. American Naturalist 155:346364.
##'
##' S.W. Kembel, M. Wu, J.A. Eisen, and J.L. Green. 2012. Incorporating 16S
##' gene copy number information improves estimates of microbial diversity and
##' abundance. PLoS Computational Biology 8(10):e1002743.
##' @keywords univar
##' @examples
##'
##' #generate random phylogeny
##' randtree <- rcoal(50)
##' #simulate trait evolution for a subset of taxa on phylogeny
##' randtraits <- sample(rTraitCont(randtree, sigma=10, root.value=100), 40)
##' #estimate trait values for "missing" taxa using PIC method
##' phyEstimate(randtree, randtraits, method="pic")
##'
##' @export phyEstimate
phyEstimate <- function(phy, trait, method = "pic", ...) {
  # trait should be a data.frame or vector with (row)names matching phylogeny
  if (is.vector(trait)) {
    trait <- data.frame(trait)
  }
  
  trait.orig <- trait
  
  # given a tree with a novel species on it
  sppObs <- row.names(trait)
  # (novel spp. are in tree but have no trait value)
  sppUnobs <- phy$tip.label[!(phy$tip.label %in% sppObs)]
  
  res <-
    as.data.frame(matrix(
      nrow = length(sppUnobs),
      ncol = 2,
      dimnames = list(sppUnobs, c("estimate", "se"))
    ))
  
  for (i in sppUnobs) {
    # for each novel species, prune all but measured + that species
    tree <- drop.tip(phy, subset(sppUnobs, sppUnobs != i))
    
    # root the tree at the novel species (leave root as trichotomy)
    tree <- root(tree, i, resolve.root = FALSE)
    
    # record branch length leading to novel species in rerooted tree
    edge <- Nnode(tree) - 1 + which(tree$tip.label == i)
    bl <- tree$edge.length[edge]
    
    # prune novel species and match new pruned tree <-> trait data
    tree <- drop.tip(tree, i)
    trait <- trait.orig[tree$tip.label,]
    
    # use PIC framework to estimate trait value at root node + error
    est <- ace(trait, tree, method = method, ...)
    val <- est$ace[1]
    cimax <- est$CI95[1, 2]
    se <- abs(cimax - val) / 1.96
    
    se.adj <- sqrt(bl) + se
    
    res[i,] <- data.frame(estimate = val, se = se.adj)
  }
  
  return(res)
}


# for discrete traits
phyEstimateDisc <-
  function(phy,
           trait,
           best.state = TRUE,
           cutoff = 0.5,
           ...) {
    # trait should be a data.frame or vector with names matching phylogeny
    
    if (is.vector(trait) || is.factor(trait)) {
      trait <- data.frame(trait)
      
      trait[, 1] <- factor(trait[, 1])
      trait.orig <- trait
      
      # given a tree with a novel taxa on it (taxa with no trait value)
      sppObs <- row.names(trait)
      sppUnobs <- phy$tip.label[!(phy$tip.label %in% sppObs)]
      trtlevels <- levels(trait[, 1])
      res <-
        as.data.frame(matrix(
          nrow = length(sppUnobs),
          ncol = length(trtlevels),
          dimnames = list(sppUnobs, trtlevels)
        ))
      
      # estimate support for different states for each novel taxon
      for (i in sppUnobs) {
        # for each novel species, prune all but measured + that species
        tree <- drop.tip(phy, subset(sppUnobs, sppUnobs != i))
        
        # root the tree at the novel species (leave root as trichotomy)
        tree <- root(tree, i, resolve.root = FALSE)
        
        
        
        
        # prune novel species and match new pruned tree <-> trait data
        tree <- drop.tip(tree, i)
        trait <- trait.orig[tree$tip.label, ]
        
        # calculate value at root node and impute to novel species
        est <- ace(trait, tree, type = "discrete", ...)
        val <- est$lik.anc[1, ]
        
        res[i, ] <- val
      }
      
      # estimate the best-supported state for each taxon
      if (best.state) {
        beststate <- as.data.frame(matrix(nrow = dim(res)[1], ncol = 2))
        colnames(beststate) <-
          c("estimated.state", "estimated.state.support")
        rownames(beststate) <- rownames(res)
        
        for (i in 1:dim(res)[1]) {
          #if >=cutoff % taxa have same label assign a consensus taxon to node
          
          best <- -sort(unlist - (res[i, ]))[1]
          print(res[i, ])
          
          ##best <- -sort(-(res[i, ]))[1]
          
          if (best >= cutoff) {
            beststate[i, 1] <- names(best)
            beststate[i, 2] <- best
          }
          else
          {
            beststate[i, 1] <- NA
            beststate[i, 2] <- NA
          }
        }
      }
      
      #return the output
      if (best.state) {
        return(cbind(as.matrix(res), beststate))
      } else {
        return(as.matrix(res))
      }
    }
  }
