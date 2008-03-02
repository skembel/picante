df2vec <- function(x, colID=1) {
	vec <- x[,colID]
	names(vec) <- row.names(x)
	vec
}
