`species.dist` <-
function (samp, metric=c("cij","jaccard","roij")) {
	metric <- match.arg(metric)
    switch(metric,
    	cij = cij(sortColumns(samp)),
    	jaccard = ( 1 - vegdist(t(sortColumns(samp)), method = "jaccard")),
    	roij = roij(sortColumns(samp))
       	)
}

