`contrast.cor.table` <-
function(nodes,cor.method="pearson") {
	concorr <- list()
	concorr$r <- cor(reflect.contrasts(nodes),method=cor.method)
	concorr$df <- length(nodes[,1])-1
	concorr$P <- contrast.cor.P(concorr$r,concorr$df)
	concorr
}

