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
	phylo.dist <- as.dist(cophenetic(sample.prune(samp,phylo)))
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
