`nri` <-
function (samp, phylo.dist, null.model = c("taxa.labels", "sample.pool", 
    "phylogeny.pool", "weighted.sample.pool"), runs = 99) 
{
    mpd.obs <- mpd(samp, phylo.dist)
    null.model <- match.arg(null.model)
    mpd.rand <- switch(null.model,
    	taxa.labels = t(replicate(runs, mpd(samp, taxaShuffle(phylo.dist)))),
    	sample.pool = t(replicate(runs, mpd(randomizeSampleKeepSampRichness(samp), phylo.dist))),
    	phylogeny.pool = t(replicate(runs, mpd(randomizeSampleKeepSampRichness(samp),
    		taxaShuffle(phylo.dist)))),
    	weighted.sample.pool = t(replicate(runs, mpd(randomizeSpeciesMatrix(samp,
    		keepSppFreq = TRUE), phylo.dist))))
    mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
    data.frame(ntaxa=specnumber(samp),mpd.obs, mpd.rand.mean, mpd.rand.sd, mpd.obs.rank, 
        mpd.obs.z, mpd.obs.p=mpd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
}

