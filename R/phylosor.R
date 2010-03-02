phylosor <- function (samp, tree) 
{
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute pd")
    }
    
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for phylosor calculation")
    }
    
    samp <- as.matrix(samp)
    s <- nrow(samp)
    phylodist <- matrix(NA, s, s)
    rownames(phylodist) <- rownames(samp)
    colnames(phylodist) <- rownames(samp)
        
    samp_comb<-matrix(NA,s*(s-1)/2,ncol(samp))
    colnames(samp_comb)<-colnames(samp)
    
    i<-1
    for (l in 1:(s - 1))
    {
    	for (k in (l + 1):s)
    	{samp_comb[i,]<-samp[l, ] + samp[k, ]
    	i<-i+1}
    }
        	
    pdsamp<-pd(samp, tree)
	pdsamp_comb<-pd(samp_comb, tree)
    
    i<-1
    for (l in 1:(s - 1)) {
        pdl <- pdsamp[l,"PD"]
        for (k in (l + 1):s) {
            pdk <- pdsamp[k,"PD"]
            pdcomb <- pdsamp_comb[i,"PD"]
            pdsharedlk <- pdl + pdk - pdcomb
            phylodist[k, l] = 2 * pdsharedlk/(pdl + pdk)
            i<-i+1
        }
    }
    return(as.dist(phylodist))
}

phylosor.rnd <- function(samp, tree, cstSor=TRUE, null.model=c("taxa.labels", "frequency", "richness", "independentswap", "trialswap"), runs=999, iterations=1000)

{
	
	Res=list()
			
	if (cstSor==TRUE)
	{
		if (null.model=="taxa.labels")
		{
			for (r in 1:runs)
			{
				Res<-c(Res,list(.phylosor.taxaShuffle(samp,tree)))}
			}
			
			else if (null.model=="richness")
			{
				for (r in 1:runs)
				{Res<-c(Res,list(.phylosor.richness(samp,tree)))}
				}
				
				else stop("This null model does not maintain Sorensen similarity: use cstSor=FALSE, or choose an other null model")
				}
	
	else
	{
		if (null.model=="taxa.labels") 
		{
			warning("This null model maintains Sorensen similarity")
			for (r in 1:runs)
			{
				Res<-c(Res,list(.phylosor.taxaShuffle(samp,tree)))
			}
			}
			
			else
			for (r in 1:runs)
			{
				Res<-c(Res,list(phylosor(randomizeMatrix(samp, null.model),tree)))
			}
			}

return(Res)
}
	

##########################################################################################
.phylosor.taxaShuffle=function(samp,tree)
	{
		sampr=samp
		colnames(sampr)=sample(colnames(samp))
		return(phylosor(sampr,tree))
		}

##########################################################################################
.phylosor.richness=function(samp,tree)
{
	s=nrow(samp)
	phylodist=matrix(NA,s,s)
	rownames(phylodist)=rownames(samp)
	colnames(phylodist)=rownames(samp)
	
	for (l in 1:(s-1))
	{
		for (k in (l+1):s)
		{
			sampr=samp
			colnames(sampr)=sample(colnames(samp))
			pdl=pd(sampr[l,,drop=FALSE],tree)$PD
			pdk=pd(sampr[k,,drop=FALSE],tree)$PD
			pdtot=pd((sampr[l,,drop=FALSE]+sampr[k,,drop=FALSE]),tree)$PD
			pdsharedlk=pdl+pdk-pdtot
			phylodist[k,l]=2*pdsharedlk/(pdl+pdk)
		}
	}
	return(as.dist(phylodist))
}			
