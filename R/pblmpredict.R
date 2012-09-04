pblmpredict<-function(x,tree1.w.novel=NULL,tree2.w.novel=NULL,predict.originals=FALSE)
{
  if (!identical(class(x),"pblm")) stop("x must be of class pblm")
  if(is.null(x$phylocovs$V1) | is.null(x$phylocovs$V2)) stop("a pblm fit with phylogenies must be supplied")
  
  sppnames1.orig<-rownames(x$phylocovs$V1)
  sppnames2.orig<-rownames(x$phylocovs$V2)
  
  if(predict.originals)
  {
    novel.spp1<-sppnames1.orig
    novel.spp2<-sppnames2.orig
    n.novels1<-length(novel.spp1)
    n.novels2<-length(novel.spp2)
    cors.1=NULL
    obs.novels1<-NULL
    pairnames<-rownames(x$variates)
    for(i in 1:n.novels1)
    {
      novel.assocs<-grep(novel.spp1[i],pairnames)
      A<-x$variates[-1*novel.assocs,]
      other.assocs<-1:length(pairnames)
      other.assocs<-other.assocs[-1*novel.assocs]
      Vyy<-V[other.assocs,other.assocs]
		  Vxy<-V[novel.assocs,other.assocs]
		  invVyy<-qr.solve(Vyy)
		  U<-rep(1,length(other.assocs))
      mu<-solve((t(U)%*%invVyy%*%U),(t(U)%*%invVyy%*%A))
		  dxgy<-Vxy%*%invVyy%*%(A-mu)
		  estmu<-rep(mu,length(dxgy))+dxgy
		  Vxgy<-V[novel.assocs,novel.assocs]-Vxy%*%invVyy%*%t(Vxy)
  
		  cors.1<-c(cors.1,cor(estmu,x$variates[novel.assocs,]))
      obs.novels1<-rbind(obs.novels1,cbind(estmu,data.frame(diag(Vxgy)),x$variates[novel.assocs,]))
    }
    names(cors.1)<-sppnames1.orig
    
  cors.2=NULL
  obs.novels2<-NULL
    for(i in 1:n.novels2)
    {
      novel.assocs<-grep(novel.spp2[i],pairnames)
      A<-x$variates[-1*novel.assocs,]
      other.assocs<-1:length(pairnames)
      other.assocs<-other.assocs[-1*novel.assocs]
      Vyy<-V[other.assocs,other.assocs]
		  Vxy<-V[novel.assocs,other.assocs]
		  invVyy<-qr.solve(Vyy)
  
  
		  U<-rep(1,length(other.assocs))
		  
      mu<-solve((t(U)%*%invVyy%*%U),(t(U)%*%invVyy%*%A))

		  dxgy<-Vxy%*%invVyy%*%(A-mu)
		  estmu<-rep(mu,length(dxgy))+dxgy

		  Vxgy<-V[novel.assocs,novel.assocs]-Vxy%*%invVyy%*%t(Vxy)
      cors.2<-c(cors.2,cor(estmu,x$variates[novel.assocs,]))
		  obs.novels2<-rbind(obs.novels2,cbind(estmu,data.frame(diag(Vxgy)),x$variates[novel.assocs,]))
    }
    names(cors.2)<-sppnames2.orig
    return(list(cors.1=cors.1[sort(names(cors.1))],cors.2=cors.2[sort(names(cors.2))],obs.novels1=obs.novels1,obs.novels2=obs.novels2))
  
  } else {
  
  A<-x$variates[,1]
  #n.cov<-dim(x$variates)[2]-1
  #full.coef<-x$coefficients[grep("full",rownames(x$coefficients)),2]
  #predicted.orig<-x$predicted
  
  if(is(tree1.w.novel)[1]=="phylo")
  {
    if(is.null(tree1.w.novel$edge.length)){tree1.w.novel<-compute.brlen(tree1.w.novel, 1)}  #If phylo has no given branch lengths
    V1<-vcv.phylo(tree1.w.novel,corr=TRUE)
  } else {
    V1<-tree1.w.novel
  }
  
  if(is(tree2.w.novel)[1]=="phylo")
  {
    if(is.null(tree2.w.novel$edge.length)){tree2.w.novel<-compute.brlen(tree2.w.novel, 1)}  #If phylo has no given branch lengths
    V2<-vcv.phylo(tree2.w.novel,corr=TRUE)
  } else {
    V2<-tree2.w.novel
  }
  
  sppnames1<-rownames(V1)
  sppnames2<-rownames(V2)
  
  V1<-as.matrix(V1)
  V2<-as.matrix(V2)
  d1<-x$signal.strength[1,2]
  d2<-x$signal.strength[2,2]
  
 	initV1<-V1
 	initV2<-V2
  nspp1<-dim(V1)[1]
  nspp2<-dim(V2)[1]
  
 	# tau = tau_i + tau_j where tau_i equals the node to tip distance
 	tau1<-matrix(diag(initV1),nspp1,nspp1) + matrix(diag(initV1),nspp1,nspp1)-2*initV1
 	tau2<-matrix(diag(initV2),nspp2,nspp2) + matrix(diag(initV2),nspp2,nspp2)-2*initV2

  V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
  V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
  V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
  V2<-V2/det(V2)^(1/nspp2)
  V<-kronecker(V2,V1)  
  invV<-qr.solve(V)
  
  #make names of species pairs
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp2))
  {
    for (u in 1:nspp1)
    {
      pairnames<-c(pairnames,paste(sppnames2[o],sppnames1[u],sep="-"))
    }
  }
  
  rownames(V)<-pairnames
  colnames(V)<-pairnames
  
  #####
  
    novel.spp1<-setdiff(sppnames1,sppnames1.orig)
    novel.spp2<-setdiff(sppnames2,sppnames2.orig)
    n.novels1<-length(novel.spp1)
    n.novels2<-length(novel.spp2)
  }
  
  if(n.novels1+n.novels2==0) stop("no novel species")
  
  other.assocs<-rownames(x$variates)
  Vyy<-V[other.assocs,other.assocs]
  invVyy<-qr.solve(Vyy)
  #U<-rep(1,length(other.assocs))
  #mu<-solve((t(U)%*%invVyy%*%U),(t(U)%*%invVyy%*%A))
	mu<-x$coefficients[1,2]
  	  
  if(n.novels1>0)
  {
    obs.novels1<-NULL
    for(i in 1:n.novels1)
    {
      novel.assocs<-grep(novel.spp1[i],pairnames)
		  Vxy<-V[novel.assocs,other.assocs]
  		dxgy<-Vxy%*%invVyy%*%(A-mu)
		  estmu<-rep(mu,length(dxgy))+dxgy
		  Vxgy<-V[novel.assocs,novel.assocs]-Vxy%*%invVyy%*%t(Vxy)
      obs.novels1<-rbind(obs.novels1,cbind(estmu,data.frame(diag(Vxgy))))
    }
  } else {obs.novels1=NULL}
  
  if(n.novels2>0)
  {
    obs.novels2<-NULL
    for(i in 1:n.novels2)
    {
      novel.assocs<-grep(novel.spp2[i],pairnames)
  	  Vxy<-V[novel.assocs,other.assocs]
  		dxgy<-Vxy%*%invVyy%*%(A-mu)
		  estmu<-rep(mu,length(dxgy))+dxgy
		  Vxgy<-V[novel.assocs,novel.assocs]-Vxy%*%invVyy%*%t(Vxy)
 	    obs.novels2<-rbind(obs.novels2,cbind(estmu,data.frame(diag(Vxgy))))
    }
  } else {obs.novels2=NULL}
  
  return(list(cors.1=NA,cors.2=NA,obs.novels1=obs.novels1,obs.novels2=obs.novels2)) 
}