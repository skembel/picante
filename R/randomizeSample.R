`randomizeSample` <-
function(samp, null.model=c("keepFreq","keepRichness","keepBoth")) {
	null.model <- match.arg(null.model)
	switch(null.model,
	keepFreq = data.frame(apply(samp,2,sample),row.names=row.names(samp)),
	keepRichness = t(data.frame(apply(samp,1,sample),row.names=colnames(samp))),
	keepBoth = matchSpeciesMatrix(samp,randomizeSpeciesMatrix(samp,keepSppFreq=TRUE)))
}

