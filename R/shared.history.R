`shared.history` <-
function(t) {
	if (is.null(t$ages)) t = node.age(t)
	terms = which(t$edge[,2]<t$Nnode+2)
	sumtips = sum(t$ages[terms])
	sumbl = sum(t$edge.length)
	shi = 1 - sumbl/sumtips
	return(shi)
	}

