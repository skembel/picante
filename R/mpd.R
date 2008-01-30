`mpd` <-
function(samp,phy.dist) {
	N <- dim(samp)[1]
	mpd <- vector()
	for (i in 1:N) {
		sppInSample <- names(samp[i,samp[i,]>0])
		sample.phy.dist <- phy.dist[sppInSample,sppInSample]
		mpd <- c(mpd,mean(sample.phy.dist[lower.tri(sample.phy.dist)]))
	}
	mpd
}

