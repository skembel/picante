`taxaShuffle` <-
function(x, strata) {
	#shuffle taxa labels on a matrix (usually a phylogenetic or trait distance matrix)
	#should die if not symmetric - can this be a dist matrix?
	#xdim <- dim(as.matrix(x))[1]
	#i <- permuted.index(xdim,strata)
	#x[i,i]
	rand.names <- resample(rownames(x))
	rownames(x) <- rand.names
	colnames(x) <- rand.names
	x
}

