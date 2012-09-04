`comm.phylo.cor` <-
function(samp,phylo,metric=c("cij","checkerboard","jaccard","doij"),
		null.model=c("sample.taxa.labels","pool.taxa.labels",
					"frequency","richness","independentswap","trialswap"),
					runs=999, ...)
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
	else for (run in 1:runs)
	{
		samp.dist <- species.dist(randomizeMatrix(samp,null.model,...),metric)
		results$random.corrs[run] <- cor(phylo.dist,samp.dist,use="pairwise")
	}
	results$obs.rank <- rank(as.vector(c(results$obs.corr,results$random.corrs)))[1]
	results$obs.rand.p <- results$obs.rank/(runs+1)
	results
}


`comm.phylo.qr` <-
function(samp,phylo,metric=c("cij","checkerboard","jaccard","doij"),
		null.model=c("sample.taxa.labels","pool.taxa.labels",
					"frequency","richness","independentswap","trialswap"),
					quant=0.75, runs=999, show.plot=FALSE, ...)
{

    if (!require(quantreg)) {
        stop("The 'quantreg' package is required to use this function.")
    }

	metric <- match.arg(metric)
	null.model <- match.arg(null.model)
	results <- list("obs.qr.intercept"=NA,"obs.qr.slope"=NA,"obs.qr.slope.p"=NA,"obs.rank"=NA,"runs"=runs,
			"random.qr.slopes"=vector(length=runs))
	phylo.dist <- as.dist(cophenetic(prune.sample(samp,phylo)))
	pool.phylo.dist <- as.dist(cophenetic(phylo))
	taxa.names <- rownames(as.matrix(phylo.dist))
	samp.dist <- as.dist(as.matrix(species.dist(samp,metric))[taxa.names,taxa.names])
    results$quantile <- quant
    qrres <- coef(rq(samp.dist~phylo.dist,tau=quant, na.action=na.omit))
    names(qrres) <- NULL
    results$obs.qr.intercept <- qrres[1] 
	results$obs.qr.slope <- qrres[2]
	
	if (show.plot) {
        plot(samp.dist~phylo.dist,xlab="Phylogenetic distance",ylab="Co-occurrence")   
    }
	
	if (null.model=="sample.taxa.labels") for (run in 1:runs)
	{
		phylo.dist <- as.dist(taxaShuffle(as.matrix(phylo.dist))[taxa.names,taxa.names])
		results$random.qr.slopes[run] <- coef(rq(samp.dist~phylo.dist,tau=quant,
		    na.action=na.omit))[2]
	}
	else if (null.model=="pool.taxa.labels") for (run in 1:runs)
	{
		phylo.dist <- as.dist(taxaShuffle(as.matrix(pool.phylo.dist))[taxa.names,taxa.names])
		results$random.qr.slopes[run] <- coef(rq(samp.dist~phylo.dist,tau=quant,
		    na.action=na.omit))[2]
	}
	else for (run in 1:runs)
	{
		samp.dist <- species.dist(randomizeMatrix(samp,null.model,...),metric)
		results$random.qr.slopes[run] <- coef(rq(samp.dist~phylo.dist,tau=quant,
		    na.action=na.omit))[2]
	}
	results$obs.rank <- rank(as.vector(c(results$obs.qr.slope,results$random.qr.slopes)))[1]
	results$obs.qr.slope.p <- results$obs.rank/(runs+1)

	if (show.plot) {
        abline(results$obs.qr.intercept,results$obs.qr.slope)   
        legend("topleft",paste("q",as.character(quant),sep=""))
    }

	results
}


`taxaShuffle` <-
function(x) {
    #TODO replace with vegan's permuted.index?
    if (!is.matrix(x)) x <- as.matrix(x)
	rand.names <- sample(rownames(x))
	rownames(x) <- rand.names
	colnames(x) <- rand.names
	return(x)
}

`tipShuffle` <- 
function(phy) {
    phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
    return(phy)
}


mpd <- function(samp, dis, abundance.weighted=FALSE) 
{
    N <- dim(samp)[1]
    mpd <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            if (abundance.weighted) {
                sample.weights <- t(as.matrix(samp[i,sppInSample,drop=FALSE])) %*% as.matrix(samp[i,sppInSample,drop=FALSE])
                mpd[i] <- weighted.mean(sample.dis,sample.weights)  

            }
            else {
                mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])         
            }
        }
        else{
            mpd[i] <- NA
        }
    }
    mpd
}


mntd <- function(samp, dis, abundance.weighted=FALSE)
{
	N <- dim(samp)[1]
	mntd <- numeric(N)
	for (i in 1:N) {
		sppInSample <- names(samp[i,samp[i,] > 0])
		if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample,sppInSample]
            diag(sample.dis) <- NA
            if (abundance.weighted)
            {
                mntds <- apply(sample.dis,2,min,na.rm=TRUE)
                sample.weights <- samp[i,sppInSample]
                mntd[i] <- weighted.mean(mntds, sample.weights)
            }
            else
            {
		        mntd[i] <- mean(apply(sample.dis,2,min,na.rm=TRUE))
		    }
		}
		else {
		    mntd[i] <- NA
		}
	}
	mntd
}


`ses.mpd` <-
function (samp, dis, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", 
    "phylogeny.pool", "independentswap", "trialswap"),
    abundance.weighted = FALSE, runs = 999, iterations = 1000) 
{
    dis <- as.matrix(dis)
    mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted)
    null.model <- match.arg(null.model)
    mpd.rand <- switch(null.model,
    	taxa.labels = t(replicate(runs, mpd(samp, taxaShuffle(dis), abundance.weighted=abundance.weighted))),
    	richness = t(replicate(runs, mpd(randomizeMatrix(samp,null.model="richness"), dis, abundance.weighted))),
    	frequency = t(replicate(runs, mpd(randomizeMatrix(samp,null.model="frequency"), dis, abundance.weighted))),    	
    	sample.pool = t(replicate(runs, mpd(randomizeMatrix(samp,null.model="richness"), dis, abundance.weighted))),
    	phylogeny.pool = t(replicate(runs, mpd(randomizeMatrix(samp,null.model="richness"),
    		taxaShuffle(dis), abundance.weighted))),
    	independentswap = t(replicate(runs, mpd(randomizeMatrix(samp,null.model="independentswap", iterations), dis, abundance.weighted))),
    	trialswap = t(replicate(runs, mpd(randomizeMatrix(samp,null.model="trialswap", iterations), dis, abundance.weighted)))
    )
    mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
    mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    mpd.obs.rank <- ifelse(is.na(mpd.rand.mean),NA,mpd.obs.rank)    
    data.frame(ntaxa=specnumber(samp),mpd.obs, mpd.rand.mean, mpd.rand.sd, mpd.obs.rank, 
        mpd.obs.z, mpd.obs.p=mpd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
}


`ses.mntd` <-
function (samp, dis, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", 
    "phylogeny.pool", "independentswap", "trialswap"),
    abundance.weighted = FALSE, runs = 999, iterations = 1000) 
{
    dis <- as.matrix(dis)
    mntd.obs <- mntd(samp, dis, abundance.weighted)
    null.model <- match.arg(null.model)
    mntd.rand <- switch(null.model,
    	taxa.labels = t(replicate(runs, mntd(samp, taxaShuffle(dis), abundance.weighted))),
    	richness = t(replicate(runs, mntd(randomizeMatrix(samp,null.model="richness"), dis, abundance.weighted))),
    	frequency = t(replicate(runs, mntd(randomizeMatrix(samp,null.model="frequency"), dis, abundance.weighted))),
    	sample.pool = t(replicate(runs, mntd(randomizeMatrix(samp,null.model="richness"), dis, abundance.weighted))),
    	phylogeny.pool = t(replicate(runs, mntd(randomizeMatrix(samp,null.model="richness"),
    		taxaShuffle(dis), abundance.weighted))),
    	independentswap = t(replicate(runs, mntd(randomizeMatrix(samp,null.model="independentswap", iterations), dis, abundance.weighted))),
    	trialswap = t(replicate(runs, mntd(randomizeMatrix(samp,null.model="trialswap", iterations), dis, abundance.weighted)))
    )
    mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
    mntd.obs.rank <- apply(X = rbind(mntd.obs, mntd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    mntd.obs.rank <- ifelse(is.na(mntd.rand.mean),NA,mntd.obs.rank)     
    data.frame(ntaxa=specnumber(samp),mntd.obs, mntd.rand.mean, mntd.rand.sd, mntd.obs.rank, 
        mntd.obs.z, mntd.obs.p=mntd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
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
    if(is.null(tree$edge.length)){tree<-compute.brlen(tree, 1)}  #If phylo has no given branch lengths
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,corr=TRUE)
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

  samp<-as.matrix(samp)
  
  if(is(tree)[1]=="phylo")
  {
    if(is.null(tree$edge.length)){tree<-compute.brlen(tree, 1)}  #If phylo has no given branch lengths
    tree<-prune.sample(samp,tree) 
    # Make sure that the species line up
    samp<-samp[,tree$tip.label, drop=FALSE]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,corr=TRUE)
  } else {
    Cmatrix<-tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    Cmatrix<-Cmatrix[species,species]
    samp<-samp[,colnames(Cmatrix),drop=FALSE]
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
    if(is.null(tree$edge.length)){tree<-compute.brlen(tree, 1)}  #If phylo has no given branch lengths
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,corr=TRUE)
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
    if(is.null(tree$edge.length)){tree<-compute.brlen(tree, 1)}  #If phylo has no given branch lengths
    tree<-prune.sample(samp,tree)
    # Make sure that the species line up
    samp<-samp[,tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix<-vcv.phylo(tree,corr=TRUE)
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

  obs.PSV<-mean(psv(samp,Cmatrix,compute.var=FALSE)[,1],na.rm=TRUE)

  # numbers of locations and species
  nlocations<-dim(samp)[1]
  nspecies<-dim(samp)[2]

  spp.PSVs<-NULL
  for (j in 1:nspecies)
  {
    spp.samp<-samp[,-j]
    spp.Cmatrix<-Cmatrix[-j,-j]
    spp.PSV<-mean(psv(spp.samp,spp.Cmatrix,compute.var=FALSE)[,1],na.rm=TRUE)
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

pd <- function (samp, tree, include.root = TRUE) 
{
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute pd")
    }
    species <- colnames(samp)
    SR <- rowSums(ifelse(samp > 0, 1, 0))
    nlocations = dim(samp)[1]
    nspecies = dim(samp)[2]
    PDs = NULL
    for (i in 1:nlocations) {
        present <- species[samp[i, ] > 0]
        treeabsent <- tree$tip.label[which(!(tree$tip.label %in% 
            present))]
        if (length(present) == 0) {
            PD <- 0
        }
        else if (length(present) == 1) {
            if (!is.rooted(tree) || !include.root) {
                warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
                PD <- NA
            }
            else {
                PD <- node.age(tree)$ages[which(tree$edge[, 2] == 
                  which(tree$tip.label == present))]
            }
        }
        else if (length(treeabsent) == 0) {
            PD <- sum(tree$edge.length)
        }
        else {
            sub.tree <- drop.tip(tree, treeabsent)
            if (include.root) {
                if (!is.rooted(tree)) {
                  stop("Rooted tree required to calculate PD with include.root=TRUE argument")
                }
                sub.tree.depth <- max(node.age(sub.tree)$ages)
                orig.tree.depth <- max(node.age(tree)$ages[which(tree$edge[,2] %in% which(tree$tip.label %in% present))])
                PD <- sum(sub.tree$edge.length) + (orig.tree.depth - 
                  sub.tree.depth)
            }
            else {
                PD <- sum(sub.tree$edge.length)
            }
        }
        PDs <- c(PDs, PD)
    }
    PDout <- data.frame(PD = PDs, SR = SR)
    rownames(PDout) <- rownames(samp)
    return(PDout)
}

`ses.pd` <-
function (samp, tree, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", 
    "phylogeny.pool", "independentswap", "trialswap"), runs = 999, iterations = 1000, ...) 
{
    pd.obs <- as.vector(pd(samp, tree, ...)$PD)
    null.model <- match.arg(null.model)
    pd.rand <- switch(null.model,
    	taxa.labels = t(replicate(runs, as.vector(pd(samp, tipShuffle(tree), ...)$PD))),
    	richness = t(replicate(runs, as.vector(pd(randomizeMatrix(samp,null.model="richness"), tree, ...)$PD))),
    	frequency = t(replicate(runs, as.vector(pd(randomizeMatrix(samp,null.model="frequency"), tree, ...)$PD))),    	
    	sample.pool = t(replicate(runs, as.vector(pd(randomizeMatrix(samp,null.model="richness"), tree, ...)$PD))),
    	phylogeny.pool = t(replicate(runs, as.vector(pd(randomizeMatrix(samp,null.model="richness"),
    		tipShuffle(tree), ...)$PD))),
    	independentswap = t(replicate(runs, as.vector(pd(randomizeMatrix(samp,null.model="independentswap", iterations), tree, ...)$PD))),
    	trialswap = t(replicate(runs, as.vector(pd(randomizeMatrix(samp,null.model="trialswap", iterations), tree, ...)$PD)))
    )
    pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
    pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
    pd.obs.z <- (pd.obs - pd.rand.mean)/pd.rand.sd
    pd.obs.rank <- apply(X = rbind(pd.obs, pd.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    pd.obs.rank <- ifelse(is.na(pd.rand.mean),NA,pd.obs.rank)    
    data.frame(ntaxa=specnumber(samp),pd.obs, pd.rand.mean, pd.rand.sd, pd.obs.rank, 
        pd.obs.z, pd.obs.p=pd.obs.rank/(runs+1),runs=runs, row.names = row.names(samp))
}