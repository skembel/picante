`ses.mnnd` <-
function (samp, dis, null.model = c("taxa.labels", "sample.pool", 
    "phylogeny.pool", "weighted.sample.pool"), runs = 99) 
{
    dis <- as.matrix(dis)
    mnnd.obs <- mnnd(samp, dis)
    null.model <- match.arg(null.model)
    mnnd.rand <- switch(null.model,
    	taxa.labels = t(replicate(runs, mnnd(samp, taxaShuffle(dis)))),
    	sample.pool = t(replicate(runs, mnnd(randomizeSample(samp,null.model="richness"), dis))),
    	phylogeny.pool = t(replicate(runs, mnnd(randomizeSample(samp,null.model="richness"),
    		taxaShuffle(dis)))),
    	weighted.sample.pool = t(replicate(runs, mnnd(randomizeSample(samp,
    		null.model = "both"), dis))))
    mnnd.obs.rank <- apply(X = rbind(mnnd.obs, mnnd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    mnnd.rand.mean <- apply(X = mnnd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    mnnd.rand.sd <- apply(X = mnnd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    mnnd.obs.z <- (mnnd.obs - mnnd.rand.mean)/mnnd.rand.sd
    data.frame(ntaxa=specnumber(samp),mnnd.obs, mnnd.rand.mean, mnnd.rand.sd, mnnd.obs.rank, 
        mnnd.obs.z, mnnd.obs.p=mnnd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
}

