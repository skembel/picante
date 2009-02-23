phylosor <- function(samp,tree)
{
  
    #If phylo has no given branch lengths
    if(is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
    }
    
    # Make sure that the species line up
    tree<-prune.sample(samp,tree)
    samp<-samp[,tree$tip.label]

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
	
.pdshort=function(comm,tree,keep.root=TRUE)
{

	nbspecies <- length(comm)
	species <- names(comm)
    absent <- species[comm==0]	#species not in sample
    present <- species[comm>0]    #species in sample
        if(length(present)==0)
    {
        #no species present
        PD<-0
    }
    else if(length(present)==1)
    {
        #one species present - PD = length from root to that tip
        if (!is.rooted(tree) || !keep.root) {
            stop("Rooted tree and keep.root=TRUE required to calculate PD of single-species communities")
        }
        
        PD <- node.age(tree)$ages[ which(tree$edge[,2] == 
                                    which(tree$tip.label==present))]
    }
    else if(length(absent)==0)
    {
        #all species present
        PD <- sum(tree$edge.length)
    }
    else
    {
        #subset of species present
        sub.tree<-drop.tip(tree,absent) 
        if (keep.root) {
            if (!is.rooted(tree) || !is.ultrametric(tree)) {
                stop("Rooted ultrametric tree required to calculate PD keeping root node")
            }
            #fixme - currently assumes ultrametric tree, could fix this
            sub.tree.depth <- max(node.age(sub.tree)$ages)
            orig.tree.depth <- max(node.age(tree)$ages)
            PD<-sum(sub.tree$edge.length) + (orig.tree.depth-sub.tree.depth)        
        }
        else {
            PD<-sum(sub.tree$edge.length)
        }
    }
    return(PD)
}


