##' Prune tree to match community data or trait data
##' 
##' Prune a phylogenetic tree to include only species present in a community
##' data set or with non-missing trait data
##' 
##' 
##' @aliases prune.sample prune.missing
##' @param phylo phylo object
##' @param samp Community data matrix
##' @param x Vector of trait data
##' @return Returns a pruned phylo object
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @keywords manip
##' @export prune.sample
`prune.sample` <-
  function(samp, phylo) {
    treeTaxa <- phylo$tip.label
    sampleTaxa <- colnames(samp)
    trimTaxa <- setdiff(treeTaxa, sampleTaxa)
    if (length(trimTaxa) > 0) drop.tip(phylo, trimTaxa) else phylo
  }

"prune.missing" <-
  function(x, phylo) {
    result <- list(NULL)
    treeTaxa <- phylo$tip.label
    traitTaxa <- names(na.omit(x[phylo$tip.label]))
    trimTaxa <- setdiff(treeTaxa, traitTaxa)
    if (length(trimTaxa) > 0) {
      result$tree <- drop.tip(phylo, trimTaxa)
    } else {
      result$tree <- phylo
    }
    result$data <- na.omit(x[phylo$tip.label])
    result
  }
