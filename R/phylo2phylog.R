`phylo2phylog` <-
function(phy, ...) {
    newick2phylog(write.tree(phy, multi.line = FALSE),...)
}

