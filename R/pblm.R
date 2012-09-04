pblm<-function(assocs,tree1=NULL,tree2=NULL,covars1=NULL,covars2=NULL,bootstrap=FALSE,nreps=10,maxit=10000,pstart=c(.5,.5)){

  # Make a vector of associations
  A<-as.matrix(as.vector(as.matrix(assocs)))
  data.vecs<-A
  
  #numbers of species and interactions
  nassocs<-length(A)
  nspp1<-dim(assocs)[1]
	nspp2<-dim(assocs)[2]
  sppnames1<-rownames(assocs)
  sppnames2<-colnames(assocs)
  #make names of species pairs
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp2))
  {
    for (u in 1:nspp1)
    {
      pairnames<-c(pairnames,paste(sppnames2[o],sppnames1[u],sep="-"))
    }
  }
    
  #Clean Covariates
  #If the covariate applies to both sets, then it should be in the matrix of the longer set 
  covnames<-NULL
  C1covs<-NULL
  if(is.null(covars1))
  {
    C1<-NULL
  } else {
    if(is.null(dim(covars1)))
    {
      C1<-matrix(covars1,nspp1,nspp2,byrow=FALSE)
      if(is.factor(covars1))
      {
        C1<-as.matrix(as.vector(C1))
        C1covs<-cbind(C1covs,C1)
        C1<-as.matrix(model.matrix(~as.factor(C1)-1)[,-1])
        colnames(C1)<-paste(rep("covar1",length(levels(covars1))-1),levels(covars1)[-1],sep="-")
      } else {
        C1<-as.matrix(as.vector(C1))
        C1covs<-cbind(C1covs,C1)
      }
      covnames<-c(covnames,"covar1")
    } else {
      C1<-NULL
      for(i in 1:dim(covars1)[2])
      {
        C1hold<-matrix(covars1[,i],nspp1,nspp2,byrow=FALSE)
        if(is.factor(covars1[,i]))
        {
          C1hold<-as.matrix(as.vector(C1hold))
          C1covs<-cbind(C1covs,C1hold)
          C1hold<-as.matrix(model.matrix(~as.factor(C1hold)-1)[,-1])
          colnames(C1hold)<-paste(rep(colnames(covars1)[i],length(levels(covars1[,i]))-1),levels(covars1[,i])[-1],sep="-")
          C1<-cbind(C1,C1hold)
        } else { 
          C1hold<-as.matrix(as.vector(C1hold))
          C1covs<-cbind(C1covs,C1hold)
          colnames(C1hold)<-colnames(covars1)[i]
          C1<-cbind(C1,C1hold)
        }
        covnames<-c(covnames,colnames(covars1)[i])
      }
    }
  
  data.vecs<-cbind(data.vecs,C1covs)
  }


  C2covs<-NULL
  if(is.null(covars2))
  {
    C2<-NULL
  } else {
    if(is.null(dim(covars2)))
    {
      C2<-matrix(covars2,nspp1,nspp2,byrow=TRUE)
      if(is.factor(covars2))
      {
        C2<-as.matrix(as.vector(C2))
        C2covs<-cbind(C2covs,C2)
        C2<-as.matrix(model.matrix(~as.factor(C2)-1)[,-1])
        colnames(C2)<-paste(rep("covar2",length(levels(covars2))-1),levels(covars2)[-1],sep="-")
      } else {
        C2<-as.matrix(as.vector(C2))
        C2covs<-cbind(C2covs,C2)
      }
      covnames<-c(covnames,"covar2")
    } else {
      C2<-NULL
      for(i in 1:dim(covars2)[2])
      {
        C2hold<-matrix(covars2[,i],nspp1,nspp2,byrow=TRUE)
        if(is.factor(covars2[,i]))
        {
          C2hold<-as.matrix(as.vector(C2hold))
          C2covs<-cbind(C2covs,C2hold)
          C2hold<-as.matrix(model.matrix(~as.factor(C2hold)-1)[,-1])
          colnames(C2hold)<-paste(rep(colnames(covars2)[i],length(levels(covars2[,i]))-1),levels(covars2[,i])[-1],sep="-")
          C2<-cbind(C2,C2hold)
        } else { 
          C2hold<-as.matrix(as.vector(C2hold))
          C2covs<-cbind(C2covs,C2hold)
          colnames(C2hold)<-colnames(covars2)[i]
          C2<-cbind(C2,C2hold)
        }
        covnames<-c(covnames,colnames(covars2)[i])
      }
    }
  
  data.vecs<-cbind(data.vecs,C2covs)
  }
  
  
# Make U, the combined matrix of covariates 
  U<-NULL
  if(is.null(C1) & is.null(C2))
  {
    U<-rep(1,length(A))
  } else {
    
    if(is.null(C1))
    {
      U<-rep(1,length(A))
    } else {
      U<-cbind(rep(1,length(A)),C1)
    }
    
    if(is.null(C2))
    {
      U<-U
    } else {
      U<-cbind(U,C2)
    }
  }

  # Begin to organize output
  if(is.null(dim(U)))
  {
    data.vecs<-data.frame(A)
    colnames(data.vecs)<-"associations"
  } else {    
    colnames(data.vecs)<-c("associations", covnames)
  }
  rownames(data.vecs)<-pairnames
  
  ######
  # Calculate Star Regression Coefficients
  #calculate for the star (assuming no phylogenetic correlation)
	astar<-solve((t(U)%*%U),(t(U)%*%A))
	MSETotal<-cov(A)
	s2aStar<-as.vector(MSETotal)*qr.solve((t(U)%*%U))
	sdaStar<-t(diag(s2aStar)^(.5))
	approxCFstar<-rbind(t(astar)-1.96%*%sdaStar, t(astar), t(astar)+1.96%*%sdaStar)
  Pstar<-U%*%astar
  Estar<-A-Pstar
  MSEStar<-cov(matrix(Estar))
  
  #######
  if(is.null(tree1) | is.null(tree2))
  {
    coefs<-approxCFstar
    rownames(coefs)<-c("lower CI 95%","estimate","upper CI 95%")
    colnames(coefs)<-paste("star",c("intercept",colnames(U)[-1]),sep="-")
    MSEs<-cbind(data.frame(MSETotal),data.frame(MSEStar))
    Pstar<-data.frame(Pstar)
    colnames(Pstar)<-"star"
    Estar<-data.frame(Estar)
    colnames(Estar)<-"star"
    output<-list(MSE=MSEs,signal.strength=NULL,coefficients=data.frame(t(coefs)),CI.boot=NULL,variates=data.frame(data.vecs),residuals=Estar,predicted=Pstar,bootvalues=NULL,Vfull=NULL)
    class(output)<-"pblm"
    return(output)
    
  } else {
  
    #tree1 is the phylogeny for the rows
    #tree2 is the phylogeny for the columns
    
    #Clean Trees
    if(is(tree1)[1]=="phylo")
    {
      if(is.null(tree1$edge.length)){tree1<-compute.brlen(tree1, 1)}  #If phylo has no given branch lengths
      V1<-vcv.phylo(tree1,corr=TRUE)
      V1<-V1[rownames(assocs),rownames(assocs)]
    } else {
      V1<-tree1[rownames(assocs),rownames(assocs)]
    }    
    
    if(is(tree2)[1]=="phylo")
    {
    if(is.null(tree2$edge.length)){tree2<-compute.brlen(tree2, 1)}  #If phylo has no given branch lengths
      V2<-vcv.phylo(tree2,corr=TRUE)
      V2<-V2[colnames(assocs),colnames(assocs)]
    } else {
      V2<-tree2[colnames(assocs),colnames(assocs)]
    }
  
    #Calculate Regression Coefficents for the base (assuming strict brownian motion evolution, ds=1)
    V1<-as.matrix(V1)
    V2<-as.matrix(V2)
    
    V1<-V1/det(V1)^(1/nspp1)   # scale covariance matrices (this reduces numerical problems caused by
  	V2<-V2/det(V2)^(1/nspp2)   # determinants going to infinity or zero)
  	V<-kronecker(V2,V1)  
    invV<-qr.solve(V)
    
    abase<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    MSEBase<-(t(A-U%*%abase)%*%invV%*%(A-U%*%abase))/(nassocs-1)  
    s2abase<-as.vector(MSEBase)*qr.solve(t(U)%*%invV%*%U)
  	sdabase<-t(diag(s2abase)^(.5))
    approxCFbase<-rbind(t(abase)-1.96%*%sdabase, t(abase), t(abase)+1.96%*%sdabase)
    Pbase<-t(t(U%*%abase)%*%invV)
    Ebase<-A-Pbase
  
    ###################
    # Full EGLS estimates of phylogenetic signal
    ##################
  	initV1<-V1
  	initV2<-V2
  
  	# tau = tau_i + tau_j where tau_i equals the node to tip distance
  	tau1<-matrix(diag(initV1),nspp1,nspp1) + matrix(diag(initV1),nspp1,nspp1)-2*initV1
  	tau2<-matrix(diag(initV2),nspp2,nspp2) + matrix(diag(initV2),nspp2,nspp2)-2*initV2
    
    # The workhorse function to estimate ds
    pegls<-function(parameters)
    {
      d1<-abs(parameters[1])
      d2<-abs(parameters[2])
  	
      V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
      V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
  
      V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
      V2<-V2/det(V2)^(1/nspp2)
      V<-kronecker(V2,V1)  
      invV<-qr.solve(V)
    
      a<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
      E<-(A-U%*%a)
      #MSE
      t(E)%*%invV%*%E/(nassocs-1)
    }
    # estimate d1 and d2 via Nelder-Mead method same as fminsearch in Matlab, by minimizing MSE
    est<-optim(pstart,pegls,control=list(maxit=maxit))        
    MSEFull<-est$value
  	d1<-abs(est$par[1])
  	d2<-abs(est$par[2])
  	
  	
    # Calculate EGLS coef w estimated ds 
  	V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
    V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
    V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
    V2<-V2/det(V2)^(1/nspp2)
    V<-kronecker(V2,V1)  
    invV<-qr.solve(V)
    aFull<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    s2aFull<-as.vector(MSEFull)*qr.solve(t(U)%*%invV%*%U)
  	sdaFull<-t(diag(s2aFull)^(.5))
  	approxCFfull<-rbind(t(aFull)-1.96%*%sdaFull, t(aFull), t(aFull)+1.96%*%sdaFull)
    Pfull<-t(t(U%*%aFull)%*%invV)
    Efull<-A-Pfull

    ########################################
    
    #organize output
    coefs<-cbind(approxCFfull,approxCFstar,approxCFbase)
    rownames(coefs)<-c("approx lower CI 95%","estimate","approx upper CI 95%")
    colnames(coefs)<-c(paste("full",c("intercept",colnames(U)[-1]),sep="-"),paste("star",c("intercept",colnames(U)[-1]),sep="-"),paste("base",c("intercept",colnames(U)[-1]),sep="-"))
    coefs<-t(coefs)
    CI.boot<-NULL
    MSEs<-cbind(data.frame(MSETotal),data.frame(MSEFull), data.frame(MSEStar), data.frame(MSEBase))
    residuals<-cbind(data.frame(Efull),data.frame(Estar),data.frame(Ebase))
    predicted<-cbind(data.frame(Pfull),data.frame(Pstar),data.frame(Pbase))
    rownames(residuals)<-pairnames
    rownames(predicted)<-pairnames
    colnames(predicted)<-c("full","star","base")
    colnames(residuals)<-c("full","star","base")
    phylocovs=list(V1=V1,V2=V2)
    
    ################
    #bootstrap CIs
    if(bootstrap)
    {
      Vtrue<-V
		  Atrue<-A
		  atrue<-aFull
		  dtrue<-c(d1,d2)
		  ehold<-eigen(Vtrue,symmetric=TRUE)
      L<-ehold$vectors[,nassocs:1]    #A or L
      G<-sort(ehold$values)      #D
      iG<-diag(G^-.5)    #iD

      # Construct Y = TT*A so that 
			# E{(Y-b)*(Y-b)'} = E{(TT*A-b)*(TT*A-b)'}
			#				  = T*V*T'
			#				  = I

      TT<-iG%*%t(L)
		  Y<-TT%*%Atrue
		  Z<-TT%*%U

		  res<-(Y-Z%*%atrue)	# residuals in orthogonalized space
		  invT<-qr.solve(TT)
			
      bootlist=NULL
      for (i in 1:nreps)
      {
        randindex<-sample(1:nassocs,replace=TRUE)	# vector of random indices
        #randindex=1:nassocs					# retain order
        YY<-Z%*%atrue+res[randindex]	# create new values of Y with random residuals
        A<-invT%*%YY	# back-transformed data
        pstart<-dtrue+c(0,.1)
        estRand<-optim(pstart,pegls,control=list(maxit=maxit))
        MSEFullrand<-estRand$value
        d1rand<-abs(estRand$par[1])
  	    d2rand<-abs(estRand$par[2])
  	
        # Calculate EGLS coef w estimated ds 
        V1<-(d1rand^tau1)*(1-d1rand^(2*initV1))/(1-d1rand^2)
        V2<-(d2rand^tau2)*(1-d2rand^(2*initV2))/(1-d2rand^2)
        V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
        V2<-V2/det(V2)^(1/nspp2)
        V<-kronecker(V2,V1)  
        invV<-qr.solve(V)
        arand<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
        
        bootlist<-rbind(bootlist,c(d1rand, d2rand, t(arand)))
      }
      nr<-dim(bootlist)[1]
      nc<-dim(bootlist)[2]
      
      #Calculate bootstrapped CIs
      alpha<-0.05  # alpha is always 0.05, but could change here
      conf<-NULL
      for(j in 1:nc)
      {
        bootconf<-quantile(bootlist[,j],probs = c(alpha/2, 1-alpha/2))
        conf<-rbind(conf,c(bootconf[1],bootconf[2]))
      }
      signal.strength<-data.frame(cbind(conf[1:2,1],dtrue,conf[1:2,2]))
      rownames(signal.strength)<-c("d1","d2")
      colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")

      #organize output
      CI.boot<-conf
      rownames(CI.boot)<-c("d1","d2","intercept",colnames(U)[-1])
      colnames(CI.boot)<-c("booted lower CI 95%","booted upper CI 95%")
      colnames(bootlist)<-c("d1","d2","intercept",colnames(U)[-1])
      output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=CI.boot,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=bootlist,phylocovs=phylocovs)
      class(output)<-"pblm"
      return(output)
    
    } else {
    ########
    # If bootstrapping not performed
    
    conf<-matrix(NA,2,2)
    signal.strength<-data.frame(cbind(conf[1,],c(d1,d2),conf[2,]))
    rownames(signal.strength)<-c("d1","d2")
    colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
    output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=NULL,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=NULL,phylocovs=phylocovs)
    class(output)<-"pblm"
    return(output)
    }
  }                                                                                                       
}