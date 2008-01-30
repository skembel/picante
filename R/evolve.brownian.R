`evolve.brownian` <-
function(phy,value=0,var=1) {
	x <- as.vector(t(evolve.phylo(phy,value,var)$tip.character))
	names(x) <- phy$tip.label
	return(x)
}
