## pic3 is a minor modification of the pic function from the ape
## library that returns the nodal values at each node of the 
## tree in a third column, along with the contrasts and the
## variance of the contrasts.

pic3 = function (x, phy, scaled = TRUE, var.contrasts = TRUE) 
{
    if (class(phy) != "phylo") 
        stop("object 'phy' is not of class \"phylo\"")
    if (is.null(phy$edge.length)) 
        stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1) 
        stop("'phy' is not rooted and fully dichotomous")
    if (length(x) != nb.tip) 
        stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x))) 
        stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")
    phy <- reorder(phy, "pruningwise")
    phenotype <- numeric(nb.tip + nb.node)
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    } else {
        if (all(names(x) %in% phy$tip.label)) 
            phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument \"x\" and the names of the tip labels did not match: the former were ignored in the analysis.")
        }
    }
    contr <- var.con <- numeric(nb.node)
    ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node), 
        as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]), 
        as.double(phy$edge.length), as.double(phenotype), as.double(contr), 
        as.double(var.con), as.integer(var.contrasts), as.integer(scaled),
        PACKAGE = "ape")
    contr <- ans[[7]]
    if (var.contrasts) {
        contr <- cbind(contr, ans[[8]])
        contr <- cbind(contr, ans[[6]][-(1:nb.tip)])
        dimnames(contr) <- list(1:nb.node + nb.tip, c("contrasts", 
            "variance","node.vals"))
    }
    else names(contr) <- 1:nb.node + nb.tip
    return(contr)
}