`nti` <-
function (samp, phylo.dist, null.model = c("taxa.labels", "sample.pool", 
    "phylogeny.pool", "weighted.sample.pool"), runs = 99) 
{
    mnnd.obs <- mnnd(samp, phylo.dist)
    null.model <- match.arg(null.model)
    mnnd.rand <- switch(null.model, taxa.labels = t(replicate(runs, 
        mnnd(samp, taxaShuffle(phylo.dist)))), sample.pool = t(replicate(runs, 
        mnnd(randomizeSampleKeepSampRichness(samp), phylo.dist))), 
        phylogeny.pool = t(replicate(runs, mnnd(randomizeSampleKeepSampRichness(samp), 
            taxaShuffle(phylo.dist)))),
        weighted.sample.pool = t(replicate(runs, 
            mnnd(randomizeSpeciesMatrix(samp, keepSppFreq = TRUE),
            	phylo.dist))))
    mnnd.obs.rank <- apply(X = rbind(mnnd.obs, mnnd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    mnnd.rand.mean <- apply(X = mnnd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    mnnd.rand.sd <- apply(X = mnnd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    mnnd.obs.z <- (mnnd.obs - mnnd.rand.mean)/mnnd.rand.sd
    data.frame(ntaxa=specnumber(samp),mnnd.obs, mnnd.rand.mean, mnnd.rand.sd, mnnd.obs.rank, 
        mnnd.obs.z, mnnd.obs.p=mnnd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
}

