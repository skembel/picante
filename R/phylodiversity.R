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

psv<-function(samp,tree,compute.var=TRUE){
  # Make samp matrix a pa matrix
  samp[samp>0]<-1
  
  flag=0
  if (is.null(dim(samp))) #if the samp matrix only has one site
  {
    samp<-rbind(samp,samp)
    flag=2
  }

  if(is(tree)[1]=="phylo")
  {
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,cor=TRUE)
  } else {
    Cmatrix<-tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    Cmatrix<-Cmatrix[species,species]
    samp<-samp[,colnames(Cmatrix)]
  }
  
    # numbers of locations and species
    SR<-rowSums(samp)
    nlocations<-dim(samp)[1]
    nspecies<-dim(samp)[2]
  
    ##################################
    # calculate observed PSVs
    #
    PSVs<-NULL
  
    for(i in 1:nlocations)
    {
      index<-seq(1:nrow(Cmatrix))[samp[i,]==1]	#species present
      n<-length(index)			#number of species present
      if(n>1)
      {
        C<-Cmatrix[index,index]	#C for individual locations
        PSV<-(n*sum(diag(as.matrix(C)))-sum(C))/(n*(n-1))
      } else {
        PSV<-NA
      }
        PSVs<-c(PSVs,PSV)
    }
    PSVout<-cbind(PSVs,SR)
  
  if(flag==2)
  {
    PSVout<-PSVout[-2,]
    return(PSVout)
  } else {
    if(compute.var==FALSE) 
    {
      return(data.frame(PSVout))
    } else {
    
      PSVvar<-NULL
      X<-Cmatrix-(sum(sum(Cmatrix-diag(nspecies))))/(nspecies*(nspecies-1));
      X<-X-diag(diag(X))
  
      SS1<-sum(X*X)/2
  
      SS2<-0;
      for(i in 1:(nspecies-1)){
        for(j in (i+1):nspecies){
          SS2<-SS2+X[i,j]*(sum(X[i,])-X[i,j])
        }
      }
      SS3<--SS1-SS2
  
      S1<-SS1*2/(nspecies*(nspecies-1))
      S2<-SS2*2/(nspecies*(nspecies-1)*(nspecies-2))
  
      if(nspecies==3){
        S3<-0
      }else{
      S3<-SS3*2/(nspecies*(nspecies-1)*(nspecies-2)*(nspecies-3))
      }
  
      for(n in 2:nspecies){
        approxi<-2/(n*(n-1))*(S1 + (n-2)*S2 + (n-2)*(n-3)*S3)
        PSVvar<-rbind(PSVvar, c(n, approxi))
      }
  
      vars<-rep(0,nlocations)
      PSVout<-cbind(PSVout,vars)
  
      for (g in 1:nlocations)
      {
        if(PSVout[g,2]>1)
        {
          PSVout[g,3]<-PSVvar[PSVout[g,2]-1,2]
        } else {
          PSVout[g,3]<-NA
        }
      }
      return(data.frame(PSVout))
    }
  }
}

psr <- function(samp,tree,compute.var=TRUE){
  PSVout<-psv(samp,tree,compute.var)
  if(is.null(dim(PSVout))==TRUE) 
  {
    PSRout<-data.frame(cbind(PSVout[1]*PSVout[2],PSVout[2]))
    names(PSRout)<-c("PSR","SR")
    rownames(PSRout)<-""
    return(PSRout)
  } else {
    PSRout<-cbind(PSVout[,1]*PSVout[,2],PSVout[,2])
    if(compute.var!=TRUE) {
      colnames(PSRout)<-c("PSR","SR")
      rownames(PSRout)<-rownames(PSVout)
      return(data.frame(PSRout))
    } else {
      PSRout<-cbind(PSRout,PSVout[,3]*(PSVout[,2])^2)       
      colnames(PSRout)<-c("PSR","SR","vars")
      rownames(PSRout)<-rownames(PSVout)
      return(data.frame(PSRout))
    }
  }
}

pse<-function(samp,tree){
  flag=0
  if (is.null(dim(samp))) #if the samp matrix only has one site
  {
    samp<-rbind(samp,samp)
    flag=2
  }

  if(is(tree)[1]=="phylo")
  {
    tree<-prune.sample(samp,tree) 
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,cor=TRUE)
  } else {
    Cmatrix<-tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    Cmatrix<-Cmatrix[species,species]
    samp<-samp[,colnames(Cmatrix)]
  }
  
  # numbers of locations and species
  SR<-apply(samp>0,1,sum)
  nlocations<-dim(samp)[1]
  nspecies<-dim(samp)[2]

  #################################
  # calculate observed phylogenetic species evenness
  PSEs<-NULL
  for(i in 1:nlocations){
    index<-seq(1,ncol(Cmatrix))[samp[i,]>0]	  # species present
    n<-length(index)                          #location species richness
    if(n>1)
    {    
    C<-Cmatrix[index,index]		                #C for individual locations
    N<-sum(samp[i,])                          #location total abundance
    M<-samp[i,samp[i,]>0]                  #species abundance column
    mbar<-mean(M)                             #mean species abundance
    PSE<-(N*t(diag(as.matrix(C)))%*%M-t(M)%*%as.matrix(C)%*%M)/(N^2-N*mbar)   #phylogenetic species evenness
    } else {PSE<-NA}
    PSEs<-c(PSEs,PSE)
  }
  PSEout=cbind(PSEs,SR)
  if (flag==2)
  {
    PSEout<-PSEout[-2,]
    return(PSEout)  
  } else {
    return(data.frame(PSEout))
  }
}



psc<-function(samp,tree){
  # Make samp matrix a pa matrix
  samp[samp>0]<-1
  flag=0
  if (is.null(dim(samp))) #if the samp matrix only has one site
  {
    samp<-rbind(samp,samp)
    flag=2
  }

  if(is(tree)[1]=="phylo")
  {
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,cor=TRUE)
  } else {
    Cmatrix<-tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    Cmatrix<-Cmatrix[species,species]
    samp<-samp[,colnames(Cmatrix)]
  }

  # numbers of locations and species
  SR<-rowSums(samp)
  nlocations<-dim(samp)[1]
  nspecies<-dim(samp)[2]

  ##################################
  # calculate observed PSCs
  #
  PSCs<-NULL

  for(i in 1:nlocations)
  {
    index<-seq(1:nrow(Cmatrix))[samp[i,]==1]	#species present
    n<-length(index)			#number of species present
    if(n>1)
    {    
      C<-Cmatrix[index,index]	#C for individual locations
      diag(C)<--1
      PSC<-sum(apply(C,1,max))/n
    } else {PSC<-NA}
      PSCs<-c(PSCs,PSC)
  }
    PSCout<-cbind(PSCs,SR)
 if (flag==2)
  {
    PSCout<-PSCout[-2,]
    return(PSCout)  
  } else {
    return(data.frame(PSCout))
  }
}

psv.spp<-function(samp,tree){
  # Make samp matrix a pa matrix
  samp[samp>0]<-1
  if (is.null(dim(samp))) #if the samp matrix only has one site
  {
    samp<-rbind(samp,samp)
  }  
  if(is(tree)[1]=="phylo")
  {
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,cor=TRUE)
  } else {
    Cmatrix<-tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    Cmatrix<-Cmatrix[species,species]
    samp<-samp[,colnames(Cmatrix)]
  }
  # reduce given Cmatrix to the species observed in samp
  samp<-samp[rowSums(samp)>1,] # prune out locations with <2 species

  #cull the species that are not found in the samp set after all communities with 1 species are removed
  preval<-colSums(samp)/sum(samp)
  indexcov<-preval>0
  Cmatrix<-Cmatrix[indexcov,indexcov]
  samp<-samp[,indexcov]

  obs.PSV<-mean(psv(samp,Cmatrix,compute.var=FALSE)[,1])

  # numbers of locations and species
  nlocations<-dim(samp)[1]
  nspecies<-dim(samp)[2]

  spp.PSVs<-NULL
  for (j in 1:nspecies)
  {
    spp.samp<-samp[,-j]
    spp.Cmatrix<-Cmatrix[-j,-j]
    spp.PSV<-mean(psv(spp.samp,spp.Cmatrix,compute.var=FALSE)[,1])
    spp.PSVs<-c(spp.PSVs,spp.PSV)
  }
  spp.PSVout<-(spp.PSVs-obs.PSV)/sum(abs(spp.PSVs-obs.PSV))
  names(spp.PSVout)<-colnames(samp)
  return(spp.PSVout)
}

psd<-function(samp,tree,compute.var=TRUE){
  if (is.null(dim(samp))) #if the samp matrix only has one site
  {
    PSDout<-data.frame(c(psv(samp,tree,compute.var)[1],psc(samp,tree)[1],psr(samp,tree,compute.var)[1],pse(samp,tree)))
    names(PSDout)<-c("PSV","PSC","PSR","PSE","SR")
    return(PSDout)
  } else {
    if (compute.var==TRUE)
    {
      PSDout<-cbind(psv(samp,tree,compute.var)[,c(1,3)],psc(samp,tree)[,1],psr(samp,tree,compute.var)[,c(1,3)],pse(samp,tree))
      colnames(PSDout)<-c("PSV","var.PSV","PSC","PSR","var.PSR","PSE","SR")
    } else {
      PSDout<-cbind(psv(samp,tree,compute.var)[,1],psc(samp,tree)[,1],psr(samp,tree,compute.var)[,1],pse(samp,tree))
      colnames(PSDout)<-c("PSV","PSC","PSR","PSE","SR")
    }
    return(data.frame(PSDout))
  }
}

pd<-function(samp,tree){
  
  # Make sure that the species line up
  tree<-prune.sample(samp,tree)
  samp<-samp[,tree$tip.label]
  species<-colnames(samp)
  
  # numbers of locations and species
  SR<-rowSums(samp)
  nlocations=dim(samp)[1]
  nspecies=dim(samp)[2]

  ##################################
  # calculate observed PDs
  #
  PDs=NULL

  for(i in 1:nlocations)
  {  
    index<-species[samp[i,]==0]	#species not in sample
    if(length(index)>=(nspecies-1))
    {
      PD<-NA
    } else {
      sub.tree<-drop.tip(tree,index)
      PD<-sum(sub.tree$edge.length)
    }
      PDs<-c(PDs,PD)
  }
  
  PDout<-data.frame(cbind(PDs,SR))
  rownames(PDout)<-rownames(samp) 
  return(PDout)
}

