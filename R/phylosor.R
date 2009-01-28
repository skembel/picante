phylosor=function(samp,tree)
{
	s=nrow(samp)
	phylodist=matrix(NA,s,s)
	rownames(phylodist)=rownames(samp)
	colnames(phylodist)=rownames(samp)
	
	for (l in 1:(s-1))
	{
		pdl=.pdshort(samp[l,],tree)
		for (k in (l+1):s)
		{
			pdk=.pdshort(samp[k,],tree)
			pdtot=.pdshort((samp[l,]+samp[k,]),tree)
			pdsharedlk=pdl+pdk-pdtot
			phylodist[k,l]=2*pdsharedlk/(pdl+pdk)
			}
			}
			return(as.dist(phylodist))
}			
phylosor.rnd=function(samp,tree,cstSor=TRUE,null.model=c("taxa.labels","frequency","richness","independentswap","trialswap"),runs=999,iterations=1000)

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
				Res<-c(Res,list(phylosor(randomizeSample(samp, null.model),tree)))
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
			pdl=.pdshort(sampr[l,],tree)
			pdk=.pdshort(sampr[k,],tree)
			pdtot=.pdshort((sampr[l,]+sampr[k,]),tree)
			pdsharedlk=pdl+pdk-pdtot
			phylodist[k,l]=2*pdsharedlk/(pdl+pdk)
			}
			}
			return(as.dist(phylodist))
}			

#############################################################################################		
	
.pdshort=function(comm,tree)
{

	nbspecies=length(comm)
	species = names(comm)
	index = species[comm == 0]
        if (length(index) >= (nbspecies - 1)) 
        {PD <- NA}
        else {
            sub.tree <- drop.tip(tree, index)
            PD <- sum(sub.tree$edge.length)}
            return(PD)}


