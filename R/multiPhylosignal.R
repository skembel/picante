'multiPhylosignal' <-
function(x,tree,...) {
	trait <- x[,1]
	names(trait) <- row.names(x)
	pruned <- pruneMissing(trait,tree)
	output <- data.frame(phylosignal(pruned$data,pruned$tree,...))
	if(length(colnames(x))>1) {
		for (i in 2:length(colnames(x))) {
			trait <- x[,i]
			names(trait) <- row.names(x)
			pruned <- pruneMissing(trait,tree)
			output <- rbind(output,phylosignal(pruned$data,pruned$tree,...))
		}
	}
	data.frame(output,row.names=colnames(x))
}