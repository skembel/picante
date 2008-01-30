`matrix2sample` <-
function(z) {
	temp <- data.frame(expand.grid(dimnames(z))[1:2], as.vector(as.matrix(z)))
	temp <- temp[(temp[, 3] > 0) & !is.na(temp[, 3]), ]
	temp <- temp[sort.list(temp[, 1]), ]
	data.frame(plot=temp[, 1], abund=temp[, 3], id=temp[, 2])
}

