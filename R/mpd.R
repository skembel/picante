`mpd` <-
function(samp,dis) {
	N <- dim(samp)[1]
	mpd <- numeric(N)
	for (i in 1:N) {
		sppInSample <- names(samp[i,samp[i,]>0])
		sample.dis <- dis[sppInSample,sppInSample]
		mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
	}
	mpd
}

