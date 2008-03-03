`mnnd` <-
function(samp,dis) {
	N <- dim(samp)[1]
	mnnd <- numeric(N)
	for (i in 1:N) {
		sppInSample <- names(samp[i,samp[i,]>0])
		sample.dis <- dis[sppInSample,sppInSample]
		diag(sample.dis) <- NA
		mnnd[i] <- mean(apply(sample.dis,2,min,na.rm=TRUE))
	}
	mnnd
}

