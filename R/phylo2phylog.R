`phylo2phylog` <-
function(phy, ...) {
    if(!require(ade4)) {stop("This function requires the ade4 package")}
    newick2phylog(write.tree(phy),...)
}

