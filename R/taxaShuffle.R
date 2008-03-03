`taxaShuffle` <-
function(x) {
    #todo replace with vegan's permuted.index?
    if (!is.matrix(x)) x <- as.matrix(x)
	rand.names <- sample(rownames(x))
	rownames(x) <- rand.names
	colnames(x) <- rand.names
	x
}

