sppregs<-function(samp,env,tree=NULL,fam="gaussian"){
  
  if(is.null(tree))
  {
    cors.phylo=NULL
  } else {                # If a tree is provided
  
    if(is(tree)[1]=="phylo")
    {
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
    cors.phylo<-Cmatrix[lower.tri(Cmatrix)] #vector of pairwise phylogenetic correlations among species
    samp<-samp[,colnames(Cmatrix)]          #only those species in the phylogeny are regressed
  }
  
  nplots<-dim(samp)[1]     # number of units in the occurence data
  nspp<-dim(samp)[2]       # number of species  in the occurence data
  sppnames<-colnames(samp) # vector of species names

  #make formula for model fit to each species

  if(is.null(dim(env)))
  {
    nenvs<-1
    envnames<-names(env)
    if (is.null(envnames)){envnames<-"env"}
    formu<-paste("y~",envnames) #Make formula
  } else {
    nenvs<-dim(env)[2]       # number of environmental variables
    envnames<-colnames(env)  # vecor of env names
    formu<-paste("y~",envnames[1]) #Make formula
    for (t in 2:nenvs)
    {
      formu<-paste(formu,envnames[t],sep="+")
    }
  }

  #Fit either the logistic or the standard regressions
  spp.resids<-NULL                           #holds each species residuals y-yhat
  spp.coef<-NULL                             #holds the coefficients of each species 
  spp.se<-NULL                               #holds the se of the coefs.
  spp.fits<-NULL                             #holds the yhats
               
  if(fam=="gaussian")
  {
      for(i in 1:nspp)
      {
        y<-samp[,i]
        mod<-glm(formu,data=cbind(y,env))
        spp.resids<-cbind(spp.resids,mod$y-mod$fitted.values)
        spp.fits<-cbind(spp.fits,mod$fitted.values)
        spp.coef<-cbind(spp.coef,summary(mod)$coef[,1])
        spp.se<-cbind(spp.se,summary(mod)$coef[,2])
      }
      cors.resid<-cor(spp.resids)[lower.tri(cor(spp.resids))]  #a vector of residual correlations among species

  } else {
  
    if (!require(brglm)) {
        stop("The 'brglm' package is required to use this function with argument fam=binomial.")
    }

    samp[samp>0]<-1             #make samp a pa matrix
    for(i in 1:nspp)
    {
      y<-samp[,i]
      mod<-brglm(formu,data=data.frame(cbind(y,env)))
      spp.resids<-cbind(spp.resids,mod$y-mod$fitted.values)
      spp.fits<-cbind(spp.fits,mod$fitted.values)
      spp.coef<-cbind(spp.coef,summary(mod)$coef[,1])
      spp.se<-cbind(spp.se,summary(mod)$coef[,2])
    }
    
    # calcualte the correlations among species residuals given that they come from a binomial process
    cor.r<-matrix(0,nrow=nspp,ncol=nspp)
    for (j in 1:dim(spp.resids)[1])
    {
      invv<-matrix((spp.fits[j,]*(1-spp.fits[j,]))^(-.5),nrow=1)
      invv<-t(invv)%*%invv
      invv[invv>(10^10)]<-(10^10)
		  r<-matrix(spp.resids[j,],nrow=1)
      cor.r=cor.r+((t(r)%*%r)*invv)
    }
    RC=cor.r/dim(spp.resids)[1]
    cors.resid<-RC[lower.tri(RC)] # Observed correlations of residuals among species
  }
  
  #add names to output
  colnames(spp.coef)<-colnames(samp)
  colnames(spp.se)<-colnames(samp)
  colnames(spp.resids)<-colnames(samp)
  colnames(spp.fits)<-colnames(samp)
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp-1))
  {
    for (u in (o+1):nspp)
    {
      pairnames<-c(pairnames,paste(sppnames[o],sppnames[u],sep="-"))
    }
  }
  cors.pa<-cor(samp)[lower.tri(cor(samp))] #obs pa correlations
  names(cors.pa)<-pairnames
  names(cors.resid)<-pairnames
  if(is.null(cors.phylo)==FALSE)
  {
    names(cors.phylo)<-pairnames
  }

  correlations<-c(cor(cors.phylo,cors.pa,use="pairwise.complete.obs"),
                  cor(cors.phylo,cors.resid,use="pairwise.complete.obs"),
                  cor(cors.phylo,cors.resid-cors.pa,use="pairwise.complete.obs"))
  names(correlations)<-c("occurence-phylo","residual-phylo","change-phylo")
  return(list(family=fam,residuals=spp.resids,coefficients=spp.coef,std.errors=spp.se,correlations=correlations,
              cors.pa=cors.pa,cors.resid=cors.resid,cors.phylo=cors.phylo))

}

sppregs.plot<-function(sppreg,rows=c(1,3),cex.mag=1,x.label="phylogenetic correlations",y.label=c("occurrence correlations w/ env","occurrence correlations wo/ env","change in correlations")){
  par(mfrow=rows,las=1,cex=cex.mag)
  plot(sppreg$cors.phylo,sppreg$cors.pa,xlab=x.label,ylab=y.label[1],main=paste("cor =",round(cor(sppreg$cors.phylo,sppreg$cors.pa,use="pairwise.complete.obs"),4)))
  abline(0,0,lty=2)
  plot(sppreg$cors.phylo,sppreg$cors.resid,xlab=x.label,ylab=y.label[2],main=paste("cor =",round(cor(sppreg$cors.phylo,sppreg$cors.resid,use="pairwise.complete.obs"),4)))
  abline(0,0,lty=2)
  plot(sppreg$cors.phylo,sppreg$cors.resid-sppreg$cors.pa,xlab=x.label,ylab=y.label[3],main=paste("cor =",round(cor(sppreg$cors.phylo,sppreg$cors.resid-sppreg$cors.pa,use="pairwise.complete.obs"),4)))
  abline(0,0,lty=2)
}