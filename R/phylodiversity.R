`comm.phylo.cor` <-
function(samp,phylo,metric=c("cij","checkerboard","jaccard","roij"),
		null.model=c("sample.taxa.labels","pool.taxa.labels",
					"frequency","richness","weighted.sample.pool"),
					runs=99)
{
	metric <- match.arg(metric)
	null.model <- match.arg(null.model)
	results <- list("obs.corr"=NA,"obs.corr.p"=NA,"obs.rank"=NA,"runs"=runs,
			"obs.rand.p"=NA,"random.corrs"=vector(length=runs))
	phylo.dist <- as.dist(cophenetic(prune.sample(samp,phylo)))
	pool.phylo.dist <- as.dist(cophenetic(phylo))
	taxa.names <- rownames(as.matrix(phylo.dist))
	samp.dist <- as.dist(as.matrix(species.dist(samp,metric))[taxa.names,taxa.names])
	results$obs.corr <- cor(phylo.dist,samp.dist,use="pairwise")
	results$obs.corr.p <- cor.test(phylo.dist,samp.dist)$p.value
	if (null.model=="sample.taxa.labels") for (run in 1:runs)
	{
		phylo.dist <- as.dist(taxaShuffle(as.matrix(phylo.dist))[taxa.names,taxa.names])
		results$random.corrs[run] <- cor(phylo.dist,samp.dist,use="pairwise")
	}
	else if (null.model=="pool.taxa.labels") for (run in 1:runs)
	{
		phylo.dist <- as.dist(taxaShuffle(as.matrix(pool.phylo.dist))[taxa.names,taxa.names])
		results$random.corrs[run] <- cor(phylo.dist,samp.dist,use="pairwise")
	}
	else if (null.model=="weighted.sample.pool") for (run in 1:runs)
	{
		samp.dist <- species.dist(randomizeSample(samp,null.model="both"),metric)
		phylo.dist <- as.dist(as.matrix(pool.phylo.dist)[rownames(as.matrix(samp.dist)),
							colnames(as.matrix(samp.dist))])	
		results$random.corrs[run] <- cor(phylo.dist,samp.dist,use="pairwise")
	}
	else for (run in 1:runs)
	{
		samp.dist <- species.dist(randomizeSample(samp,null.model),metric)
		results$random.corrs[run] <- cor(phylo.dist,samp.dist,use="pairwise")
	}
	results$obs.rank <- rank(as.vector(c(results$obs.corr,results$random.corrs)))[1]
	results$obs.rand.p <- results$obs.rank/(runs+1)
	results
}


mpd <- function(samp, dis) 
{
    N <- dim(samp)[1]
    mpd <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
        }
        else{
            mpd[i] <- 0
        }
    }
    mpd
}


`mnnd` <-
function(samp,dis) {
	N <- dim(samp)[1]
	mnnd <- numeric(N)
	for (i in 1:N) {
		sppInSample <- names(samp[i,samp[i,]>0])
		if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample,sppInSample]
            diag(sample.dis) <- NA
		    mnnd[i] <- mean(apply(sample.dis,2,min,na.rm=TRUE))
		}
		else {
		    mnnd[i] <- 0
		}
	}
	mnnd
}

`ses.mpd` <-
function (samp, dis, null.model = c("taxa.labels", "sample.pool", 
    "phylogeny.pool", "weighted.sample.pool"), runs = 99) 
{
    dis <- as.matrix(dis)
    mpd.obs <- mpd(samp, dis)
    null.model <- match.arg(null.model)
    mpd.rand <- switch(null.model,
    	taxa.labels = t(replicate(runs, mpd(samp, taxaShuffle(dis)))),
    	sample.pool = t(replicate(runs, mpd(randomizeSample(samp,null.model="richness"), dis))),
    	phylogeny.pool = t(replicate(runs, mpd(randomizeSample(samp,null.model="richness"),
    		taxaShuffle(dis)))),
    	weighted.sample.pool = t(replicate(runs, mpd(randomizeSample(samp,
    		null.model = "both"), dis))))
    mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
    data.frame(ntaxa=specnumber(samp),mpd.obs, mpd.rand.mean, mpd.rand.sd, mpd.obs.rank, 
        mpd.obs.z, mpd.obs.p=mpd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
}

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
