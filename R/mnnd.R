`mnnd` <-
function(samp,phy.dist) {
	N <- dim(samp)[1]
	mnnd <- vector()
	for (i in 1:N) {
		sppInSample <- names(samp[i,samp[i,]>0])
		sample.phy.dist <- phy.dist[sppInSample,sppInSample]
		diag(sample.phy.dist) <- NA
		mnnd <- c(mnnd,mean(sapply(data.frame(sample.phy.dist),min,na.rm=TRUE)))
	}
	mnnd
}

