
###################################################################################################################
#PGLMM code from Ives A.R. & Helmus M.R. (2011). Generalized linear mixed models for phylogenetic analyses of community structure. Ecological Monographs, 81, 511-525.
#
# This code contains the functions pglmm.sim, pglmm.data,
###################################################################################################################



#############################################################################################
# Simulates models used to test the PGLMM models
# 23 June 2010 MATLAB CODE GLMM_model_simulations.m
# 6 September 2011 R CODE TRANSLATION - HELMUS
#############################################################################################

pglmm.sim <- function(tree,nsites=30,modelflag=1,figs=TRUE,second.env=TRUE,compscale = 1)
{


  bspp2 <- NULL
  Vcomp <- NULL
  envirogradflag2 <- 0

  if (is(tree)[1] == "phylo")
  {
    if (is.null(tree$edge.length))
    {
      tree <- compute.brlen(tree, 1)
    }
    V <- vcv.phylo(tree, corr = TRUE)
  } else {
    V <- tree
  }
   nspp<-dim(V)[1]
  # parameter values for each model
  if(modelflag==1 | modelflag==2)
  {
    Xscale <- 1
    Mscale <- .5
    Vscale1 <- 1
    Vscale2 <- 1
    b0scale <- .5

    envirogradflag1 <- 1
    if(second.env){envirogradflag2 <- 1} #envirogradflag2 <- 1
    elimsitesflag <- 1
    repulseflag <- 0
  }

  if(modelflag==3)
  {
    Xscale <- 1
    Mscale <- 0
    Vscale1 <- 1
    Vscale2 <- 1
    compscale <- compscale
    b0scale <- 0

    # repulsion matrix
    Vcomp <- solve(V,diag(nspp))
    Vcomp <- Vcomp/max(Vcomp)
    Vcomp <- compscale*Vcomp
    iDcomp <- t(chol(Vcomp))
    colnames(Vcomp)<-rownames(Vcomp)

    envirogradflag1 <- 1
    if(second.env){envirogradflag2 <- 1} # envirogradflag2=0;
    elimsitesflag <- 0
    repulseflag <- 1
  }

  if(modelflag==4)
  {
    Xscale <- 1
    Mscale <- .5

    Vscale1 <- .05
    Vscale2 <- 1

    b0scale <- .5

    envirogradflag1 <- 1
    if(second.env){envirogradflag2 <- 1} # envirogradflag2=0;
    elimsitesflag <- 1
    repulseflag <- 0
  }

  if(modelflag==5)
  {
    Xscale <- 1
    Mscale <- .5

    Vscale1 <- 1
    Vscale2 <- 1

    b0scale <- .5

    envirogradflag1 <- 1
    if(second.env){envirogradflag2 <- 1} #envirogradflag2=1

    elimsitesflag <- 1
    repulseflag <- 0
  }

  Vforsim <- V
  iD <- t(chol(Vforsim))

  # set up environmental gradient among sites

  mx <- t(as.matrix((-(nsites)/2):(nsites/2)))
  # number of sites
  m <- length(mx)


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #% set up independent variables
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  #% environmental gradient 1
  if(envirogradflag1 == 1)
  {
    e <- iD %*% rnorm(nspp)
    e <- Vscale1*(e-mean(e))/apply(e,2,sd)
    bspp1 <- e

    X <- 1/(1+exp(-(b0scale*array(1,c(m,1)) %*% rnorm(nspp) + t(mx) %*% t(e))))
    X <- Xscale*X
    Xsmooth <- X

    X1 <- diag(1-Mscale*runif(m)) %*% X
  }

  # environmental gradient 2
  if(envirogradflag2 == 1)
  {
    e <- iD %*% rnorm(nspp)
    e <- Vscale2*(e-mean(e))/apply(e,2,sd)
    bspp2 <- e

    mx2 <- as.matrix(mx[sample(m)])
    X <- 1/(1+exp(-(b0scale*array(1,c(m,1)) %*% rnorm(nspp) + mx2 %*% t(e))))
    X <- Xscale*X
    Xsmooth <- Xsmooth*X

    X2 <- diag(1-Mscale*runif(m)) %*% X
    X <- X1*X2
  } else {
    X <- X1
  }

  # phylogenetic repulsion
  if(repulseflag == 1)
  {
    bcomp <- NULL
    for(i in 1:m)
    {
        bcomp <- cbind(bcomp, iDcomp %*% rnorm(nspp))
    }
    bcomp0 <- 0
    Xcomp <- exp(bcomp0+bcomp)/(1+exp(bcomp0+bcomp))
    X <- X*t(Xcomp)
  }

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #% simulate data
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  #% convert distribution to presence/absence
  Y <- matrix(0,ncol=nspp,nrow=m)
  Y[matrix(runif(nspp*m),ncol=nspp) < X] <- 1
  colnames(Y)<-colnames(X)

  # eliminate sites with 0 spp
  pick <- (rowSums(Y)>0)
  Y <- Y[pick,]
  mx <- mx[pick]
  m <- length(mx)

  # eliminate sites with 1 spp
  if(elimsitesflag == 1){
    pick <- (rowSums(Y)>1)
    Y <- Y[pick,]
    mx <- mx[pick]
    m <- length(mx)
  }

  if(figs)
  {
    if (!require(plotrix))
    {
      stop("The 'plotrix' package is required to plot images from this function")
    }

    if(.Platform$OS.type == "unix") quartz() else windows()
    par(mfrow=c(5,1),las=1,mar=c(2, 4, 2, 2) - 0.1)
    matplot(Xsmooth,type="l",ylab="probability",main="Xsmooth")
    matplot(X,type="l",ylab="probability",main="X")
    hist(colSums(Y),main="spp per site")
    hist(rowSums(Y),main="sites per spp")
    plot(x=mx,y=rowSums(Y),main="SR across gradient",type="l")

    if(.Platform$OS.type == "unix") quartz() else windows()
    color2D.matplot(1-X/max(X),xlab="species",ylab="sites",main="probabilities")

    if(.Platform$OS.type == "unix") quartz() else windows()
    color2D.matplot(1-Y,xlab="species",ylab="sites",main="presence-absence")
  }

  return(list(Vphylo=V,Vcomp=Vcomp,Y=Y,X=X,u=mx,bspp1=bspp1,bspp2=bspp2))
}

######################################################################################################################
##Function that organizes the data so that PGLMM can be fit.
######################################################################################################################

pglmm.data<-function(modelflag=1,sim.dat=NULL,samp=NULL,tree=NULL,traits=NULL,env=NULL,Vcomp=NULL)
{

  if(!is.null(sim.dat))
  {
    tree<-sim.dat$Vphylo
    Vcomp<-sim.dat$Vcomp
    samp<-sim.dat$Y
    traits<-sim.dat$bspp1
    env<-sim.dat$u
  }

  is.empty<-function(x){length(x) == 0}
  if(is.empty(samp))
  {
    stop("sample matrix (Y) is empty")
  }

  if (is(tree)[1] == "phylo")
  {
    if (is.null(tree$edge.length))
    {
      tree <- compute.brlen(tree, 1)
    }
    tree <- prune.sample(samp, tree)
    samp <- samp[, tree$tip.label]
    V <- vcv.phylo(tree, corr = TRUE)
    species <- colnames(samp)
    preval <- colSums(samp)/sum(sum(samp))
    species <- species[preval > 0]
    V <- V[species, species]
    Vcomp <- Vcomp[species, species]
    samp <- samp[, colnames(V)]
    traits <- as.matrix(traits[species,])


  } else {
    V <- tree
    species <- colnames(samp)
    preval <- colSums(samp)/sum(sum(samp))
    species <- species[preval > 0]
    V <- V[species, species]
    Vcomp <- Vcomp[species, species]
    samp <- samp[, colnames(V)]
    traits <- as.matrix(traits[species,])
    }

  # X = species independent variables (traits)
  # U = site independent variables (environment)
  # Y = site (rows) by species (columns) presence/absence matrix, the dependent variable (binary 0,1)
  # V = phylogenetic covariance matrix
  Y<-samp
  U<-matrix(env,nrow=length(env),ncol=1)
  X<-traits

  nsites<-dim(Y)[1]
  nspp<-dim(Y)[2]

  # create dependent variable vector
  YY<-t(Y)
  YY<-as.matrix(as.vector(as.matrix(YY)))

  if(modelflag==1)
  {
    # set up covariance matrices
    Vfullspp<-kronecker(diag(nsites),V)   #should be nspp*nsites in dimension
    Vfullsite<-kronecker(diag(nsites),matrix(1,nspp,nspp))

    #VV<-abind(Vfullspp,Vfullsite,along=3)  #could potentially also use the abind function in the abind library but I will see if we can get this to work instead
    VV<-list(Vfullspp=Vfullspp,Vfullsite=Vfullsite)

    # create independent variables
    XX<-kronecker(matrix(1,nsites,1),diag(nspp))
    return(list(YY=YY,VV=VV,XX=XX))
  }

  if(modelflag==2)
  {
    # set up covariance matrices
    u<-scale(U)
    U<-kronecker(u,matrix(1,nspp,1))

    Vfullspp<-kronecker(matrix(1,nsites,nsites),diag(nspp))
    VfullsppV<-kronecker(matrix(1,nsites,nsites),V)

    VfullUCU<-diag(as.vector(U)) %*% Vfullspp %*% diag(as.vector(U))
    VfullUCUV<-diag(as.vector(U)) %*% VfullsppV %*% diag(as.vector(U))

    Vfullsite<-kronecker(diag(nsites),matrix(1,nspp,nspp))

	  VV<-list(VfullUCU=VfullUCU,VfullUCUV=VfullUCUV,Vfullsite=Vfullsite)

    # create independent variables
    XXspp<-kronecker(matrix(1,nsites,1),diag(nspp))
    XX<-cbind(U,XXspp)

    # create dependent variable vector
    YY<-as.vector(t(Y))

    return(list(YY=YY,VV=VV,XX=XX))
  }

  if(modelflag == 3)
  {
    # set up covariance matrices
    u<-scale(U)
    U<-kronecker(u,matrix(1,nspp,1))

    #Vfullspp<-kronecker(diag(nsites),V)

    #VfullUCU<-diag(as.vector(U)) %*% Vfullspp %*% diag(as.vector(U))

    Vfullsite<-kronecker(diag(nsites),matrix(1,nspp,nspp))
    if(is.null(Vcomp))
    {
      compscale<-1
      # repulsion matrix
      Vcomp <- solve(V,diag(nspp))
      Vcomp <- Vcomp/max(Vcomp)
      Vcomp <- compscale*Vcomp
    }
    Vfullcomp<-kronecker(diag(nsites),Vcomp)

	  VV<-list(Vfullcomp=Vfullcomp,Vfullsite=Vfullsite)

    # create independent variables
    XXspp<-kronecker(matrix(1,nsites,1),diag(nspp))
    XX<-cbind(XXspp, ((U %*% matrix(1,1,nspp))*XXspp))

    # create dependent variable vector
    YY<-as.vector(t(Y))

    return(list(YY=YY,VV=VV,XX=XX))
  }

  if(modelflag == 4)
  {
    # set up covariance matrix for traits
    Vfulltrait <- kronecker(diag(nsites),traits %*% t(traits))
    traitscale4 <- 100
    Vfulltrait <- traitscale4*Vfulltrait
    # set up covariance matrices
    Vfullsite<-kronecker(diag(nsites),matrix(1,nspp,nspp))
 	  VV<-list(Vfulltrait=Vfulltrait,Vfullsite=Vfullsite)

    # create independent variables
    XXspp<-kronecker(matrix(1,nsites,1),diag(nspp))
    XX<-XXspp

    # create dependent variable vector
    YY<-as.vector(t(Y))
    return(list(YY=YY,VV=VV,XX=XX))
  }

  if(modelflag == 5)
  {
    # set up covariance matrices
    Vfulltrait <- kronecker(diag(nsites),traits %*% t(traits))
    traitscale5 <- 10
    Vfulltrait <- traitscale5*Vfulltrait

    Vfullspp<-kronecker(diag(nsites),V)   #should be nspp*nsites in dimension

    Vfullsite<-kronecker(diag(nsites),matrix(1,nspp,nspp))

    VV<-list(Vfulltrait=Vfulltrait,Vfullspp=Vfullspp,Vfullsite=Vfullsite)

    # create independent variables
    XXspp<-kronecker(matrix(1,nsites,1),diag(nspp))
    XX<-XXspp

    # create dependent variable vector
    YY<-as.vector(t(Y))
    return(list(YY=YY,VV=VV,XX=XX))
  }
}

###########################################################################################################
# Calls a PGLMM for logistic regression, workhorse function
###########################################################################################################
#global tH tinvW tVV tX
tH<<-NULL
tinvW<<-NULL
tVV<<-NULL
tX<<-NULL
pglmm.fit<-function(dat=NULL,Y=NULL,X=NULL,VV=NULL,sp.init=0.5,maxit=25,exitcountermax=50) # [B,B0,s,S95int,LL,flag]=PGLMM_nosparse_funct(Y,X,VV,s)
{
  #####################################################################################################
  #########REML estimation Function not called directly by the user
  #####################################################################################################
  PGLMM.reml<-function(sp)
  {


    sp<-abs(Re(sp))
    Cd<-matrix(0,dim(tVV[[1]])[1],dim(tVV[[1]])[2])
    for(i in 1:length(sp))
    {
	   Cd=Cd + sp[i] * tVV[[i]]
    }
    V<- tinvW + Cd
    invV<- solve(V,diag(x=1,nrow=dim(V)[1],ncol=dim(V)[2]))

    # check to see if V is positive definite
    if(all( eigen(V)$values >0 ))
    {
      cholV<-chol(V)
      # ML
      # LL=.5*(2*sum(log(diag(chol(V))))+tH'*invV*tH)

      #REML
      LL=.5 * (2 * sum(log(diag(cholV))) + t(tH) %*% invV %*% tH + log(det(t(tX) %*% invV %*% tX)))
    } else {
      LL<-10^10
    }
    LL
  }


    if (!require(corpcor))
    {
      stop("The 'corpcor' package is required")
    }

  # dat = list from PGLMM.data
  # X = independent variable
  # Y = dependent variable (binary 0,1)
  # VV = list containing covariance matrices
  # s = initial values of s

  # Returns
  # B_SE = coefficients with SEs from GLMM
  # B0_SE = coefficients with SEs from logistic regression
  # s = parameter(s) of the covariance matrix
  # flag = 1 if convergence achieved, 0 otherwise

  #global tH tinvW tVV tX Lmin S iss
  if(!is.null(dat))
  {
    X<-dat$XX
    Y<-dat$YY
    VV<-dat$VV
  }

  is.empty<-function(x){length(x) == 0}
  if(any(unlist(lapply(list(X=X,Y=Y,VV=VV),is.empty))))
  {
    stop("a data matrix is empty")
  }

  if(any(unlist(lapply(list(X=X,Y=Y,VV=VV),is.na))))
  {
    stop("a data matrix is NA")
  }

  n<-dim(X)[1]
  p<-dim(X)[2]

  #initial estimates for the s parameters (scaling parameters of the covariance matrices
  sp<-matrix(sp.init,length(as.list(VV)),1)

  B0<-matrix(mean(Y),p,1)
  oldB0<-matrix(10,p,1)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #% initial values for B
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  counter<-0
  #sparse<-function(W){as(Matrix(W),"sparseMatrix")}
  while (matrix((t(oldB0-B0) %*% (oldB0-B0) > 10^-8)) & (counter<100))
  {
    oldB0<-B0
    counter<-counter+1

    mu<-exp(X %*% B0) / (1+exp(X %*% B0))
    W<-as.vector((mu*(1-mu))^-1)
    invW<-(diag(W))
    #invW<-sparse(diag(W))

    Z<-(X %*% B0 + (Y-mu)/(mu*(1-mu)))
    denom<-(t(X) %*% invW %*% X)
    #Z<-sparse(X %*% B0 + (Y-mu)/(mu*(1-mu)))
    #denom<-sparse(t(X) %*% invW %*% X)

    options(warn=-1)
    if(any(c(is.nan(denom),is.infinite(denom),is.null(denom))))
    {
		  B0<-solve((t(X)%*% X),(t(X) %*% Y))
      (counter<-100)
    } else {
      #num<-sparse(t(X) %*% invW %*% Z)
      num<-(t(X) %*% invW %*% Z)
      B0<-pseudoinverse(denom) %*% num
    }
  }

  if (is.nan(B0) | is.infinite(B0))
  {
	 mu<-matrix(mean(Y),p,1)
	 B0<-log(mu/(1-mu))
  }

  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## GLMM
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # initial estimates
  B<-B0
  b<-matrix(0,n,1)
  #bet<-sparse(rBind(B,b))
  #mu<-sparse(exp(X %*% B)/(1+exp(X %*% B)))
  bet<-(rbind(B,b))
  mu<-(exp(X %*% B)/(1+exp(X %*% B)))  # the mean function for the binomial process change this to get at other distributions

  # matrix including fixed covariates and dummies for "random effects"
  XX<-cbind(X,diag(n))
  Cdum<-matrix(0,n,n)
  for(i in 1:length(sp))
  {
	 Cdum=Cdum+(sp[i]*VV[[i]])
  }
  #Cdum<-Matrix(Cdum)
  est<-t(rbind(sp,B))
  oldest<-matrix(10^6,dim(est)[1],dim(est)[2])

  exitcounter<-0
  while (matrix(((est-oldest) %*% t(est-oldest)) > 10^-4) & exitcounter<=exitcountermax)
  {
    oldest<-est
    # Using notation of Breslow & Clayton (1993)
    # Note: Schall (1991) equation for V incomplete; see McCullagh & Nelder p. 40 and
    # Breslow & Clayton (1993)
    W<-as.vector((mu*(1-mu))^-1)
    #invW<-sparse(diag(W))
    invW<-(diag(W))

    V<-invW+Cdum  #Breslow & Clayton (1993) between eq 10 and 11
    invV<-solve(V,diag(n))  # needed for eq 10

    #################### I did not deal with this code for when the solve function returns a singular matrix  ###################
    #ww=lastwarn;
    #lastwarn('noper')
    #if sum(ww(1:5)=='Matri')==5
		#'Matrix close to singular: alternative B0 tried'
		#mu=rand(length(X(1,:)),1);
		#B0=log(mu./(1-mu));
		#b=zeros(n,1);
		#beta=[B0;b];

		#mu=exp(X*B0)./(1+exp(X*B0));

		#invW=diag((mu.*(1-mu)).^-1);
		#V=invW+C;
		#invV=V\eye(n);

		#B0'
    #end
    #########################################

    Z<-X %*% B + b + (Y-mu)/(mu*(1-mu))
    denom<-(t(X) %*% invV %*% X)    # left side of eq 10
    #denom<-sparse(t(X) %*% invW %*% X)
    num<-(t(X) %*% invV %*% Z) #right side of eq 10
    B<-pseudoinverse(denom) %*% num  # solve eq 10 for the fixed effects (alpha in Breslow & Clayton (1993))
    b<-Cdum %*% invV %*% (Z-X %*% B) #eq 11

    bet<-(rbind(B,b))
    mu<-exp(XX %*% bet)/(1+exp(XX %*% bet))

    #DEFINE THESE AS GLOBAL tH tinvW tVV tX
    tH<-Z - X %*% B
    tinvW<-diag(as.vector((mu*(1-mu))^-1))
    tX<-X
    tVV<-VV

    # call to obtain estimates of covariance matrix parameters s
    #options=optimset('MaxFunEvals',25,'MaxIter',25,'TolX',10^-2,'TolFun',10^-2);
    #pars<-list(sp=sp,tVV=tVV,tH=tH,tinvW=tinvW,tVV=tVV,tX=tX)
    sp<-abs(optim(sp,PGLMM.reml,control=list(maxit=maxit,abstol=10^-1))$par)

    if(exitcounter==0){
      scrn.output<-c(exitcounter,t(sp), t(B[1:4]))
      names(scrn.output)<-c("exitcounter",paste("sigma",1:length(sp)),paste("B",1:4))
      print(scrn.output)
    } else {
      print(c(exitcounter,t(sp), t(B[1:4])))
    }

    Cdum<-matrix(0,n,n)
    for(i in 1:length(sp))
    {
	     Cdum=Cdum + sp[i] * tVV[[i]]
    }
    est<-t(rbind(sp,B))
    exitcounter<-exitcounter+1
  }

  # flag cases of non-convergence
  flag <- "converged"
  if(exitcounter>=exitcountermax)
  {
	 flag<-"did not converge, try increasing exitcountermax"
  }

  ## flag cases of non-convergence
  #if(is.nan(B)){
	# return
  #}

  ##############
  # compute final estimates and SEs
  W<-as.vector((mu*(1-mu))^-1)
  #invW<-sparse(diag(W))
  invW<-(diag(W))
  V<-invW+Cdum
  invV<-solve(V,diag(n))

  Z<-X %*% B + b + (Y-mu)/(mu*(1-mu))
  B<-solve((t(X) %*% invV %*% X),(t(X) %*% invV %*% Z))

  Chi025<-5.0238
  Chi05<-3.8414

  #options=optimset('MaxFunEvals',100,'MaxIter',100,'TolX',5*10^-3,'TolFun',10^-4);

  S95int<-NULL
  S<-sp
  for(iss in 1:length(sp))
  {
    Lmin<-PGLMM.reml(S)
    #Smin
    if(sp[iss] > 0.02)
    {Smin<-optimize(PGLMM.reml,c(0,.9*sp[iss]),tol = 0.0001)$minimum
    } else {
    Smin<-0}
    if(Smin<0.01){Smin<-0}
    #Smax
	   if(sp[iss] > 0.02){
	     Smax<-optimize(PGLMM.reml,c(1.1*sp[iss],10*sp[iss]),tol = 0.0001)$minimum
	   } else {
	     Smax<-optimize(PGLMM.reml,c(.1,5),tol = 0.0001)$minimum}
    S95int<-rbind(S95int,c(abs(Smin),abs(Smax)))
  }

  names(sp)<-names(VV)
  colnames(S95int)<-c("0.05","0.95")
  return(list(B=B,B0=B0,s=cbind(sp,S95int),LL=Lmin,flag=flag))
}


###########################END