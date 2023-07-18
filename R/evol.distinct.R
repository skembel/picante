# Evolutionary distinctiveness by:
# a) equal splits (Redding and Mooers 2006)
# b) fair proportions (Isaac et al., 2007)
# The scale option refers to whether or not the phylogeny should be scaled
# to a depth of 1 or, in the case of an ultrametric tree,  scaled such that
# branch lengths are relative.
# If use.branch.lengths=FALSE, then all branch lengths are changed to 1.



##' Species' evolutionary distinctiveness
##' 
##' Calculates evolutionary distinctiveness measures for a suite of species by:
##' a) equal splits (Redding and Mooers 2006) b) fair proportions (Isaac et
##' al., 2007). Returns a datafram with species identifiers and species scores.
##'
##' @examples
##'
##'data(phylocom)
##'evol.distinct(phylocom$phylo, type = 'equal.splits')
##'
##' 
##' @param tree an object of class phylo
##' @param type a) equal splits (Redding and Mooers 2006) or b) fair
##' proportions (Isaac et al., 2007)
##' @param scale The scale option refers to whether or not the phylogeny should
##' be scaled to a depth of 1 or, in the case of an ultrametric tree, scaled
##' such that branch lengths are relative.
##' @param use.branch.lengths If use.branch.lengths=FALSE, then all branch
##' lengths are changed to 1.
##' @note This function will return a vector of evolutionary distinctivenss for
##' every species in the given tree. If only a subset of values are needed
##' there are two, concetually distinct options: either prune the tree first
##' and then pass the tree in or subset the resulting vector.  These two
##' options will provide very different outputs.
##' @author Karen Magnuson-Ford, Will Cornwell, Arne Mooers, Mark Vellend
##' @references Redding, D.W. and Mooers, A.O. (2006). Incorporating
##' evolutionary measures into conservation prioritisation. Conservation
##' Biology, 20, 1670-1678.
##' 
##' Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M.
##' (2007). Mammals on the EDGE: conservation priorities based on threat and
##' phylogeny. PLoS ONE, 2, e296.
##' 
##' Mark Vellend, William K. Cornwell, Karen Magnuson-Ford, and Arne Mooers. In
##' press. Measuring phylogenetic biodiversity. In: Biological diversity:
##' frontiers in measurement and assessment.  Edited by Anne Magurran and Brian
##' McGill.
##' @export evol.distinct
evol.distinct <- function(tree, type = c("equal.splits", "fair.proportion"),
                          scale = FALSE, use.branch.lengths = TRUE) {
  type <- match.arg(type)

  if (is.rooted(tree) == FALSE) {
    warning("A rooted phylogeny is required for meaningful
output of this function", call. = FALSE)
  }

  if (scale == TRUE) {
    # Scale tree to have unit depth (for an ultrametric tree)
    # or scale all branches to unit length (for an additive tree)

    if (is.ultrametric(tree) == TRUE) {
      tree$edge.length <- tree$edge.length /
        (as.numeric(branching.times(tree)[1]))
    } else {
      tree$edge.length <- tree$edge.length / sum(tree$edge.length)
    }
  }

  if (use.branch.lengths == FALSE) {
    tree$edge.length <- rep(1, length(tree$edge.length))
  }


  for (i in seq_along(tree$tip.label)) {
    spp <- tree$tip.label[i]
    nodes <- .get.nodes(tree, spp)
    # get rid of root node
    nodes <- nodes[1:(length(nodes) - 1)]

    internal.brlen <- tree$edge.length[which(tree$edge[, 2] %in% nodes)]

    # apportion internal branch lengths appropriately
    if (length(internal.brlen) != 0) {
      internal.brlen <- internal.brlen *
        switch(type,
          equal.splits = sort(rep(
            0.5,
            length(internal.brlen)
          )^seq_along(internal.brlen)),
          fair.proportion = {
            for (j in seq_along(nodes)) {
              sons <- .node.desc(tree, nodes[j])
              n.descendents <- length(sons$tips)
              if (j == 1) {
                portion <- n.descendents
              } else {
                portion <- c(n.descendents, portion)
              }
            }
            1 / portion
          }
        )
    }

    # sum internal branch lengths with the pendant edge
    ED <- sum(internal.brlen, tree$edge.length[which.edge(tree, spp)])

    if (i == 1) {
      w <- ED
    } else {
      w <- c(w, ED)
    }
  }
  results <- cbind(tree$tip.label, as.data.frame(w))
  names(results) <- c("Species", "w")
  return(results)
}
