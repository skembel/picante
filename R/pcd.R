#Phylogenetic Community Dissimilarity from
#Ives A.R. & Helmus M.R. (2010). Phylogenetic metrics of community similarity. The American Naturalist, 176, E128-E142.

pcd <- function(comm, tree, PSVmncd=NULL, PSVpool=NULL, reps=10^4)
{
  SSii<-PSVmncd
  SCii<-PSVpool

  # Make comm matrix a pa matrix
  comm[comm>0]<-1

	# convert trees to VCV format
	if (is(tree)[1] == "phylo")
  {
		if (is.null(tree$edge.length)) {tree <- compute.brlen(tree, 1)}    #If phylo has no given branch lengths
		tree <- prune.sample(comm, tree)
		V <- vcv.phylo(tree, corr = TRUE)
		comm <- comm[, tree$tip.label]
	} else {
		V <- tree
    species<-colnames(comm)
    preval<-colSums(comm)/sum(comm)
    species<-species[preval>0]
    V<-V[species,species]
    comm<-comm[,colnames(V)]
  }

  if (!is.null(SSii) & length(SSii)!=max(rowSums(comm))) {
        stop("The length of PSVmncd is not equal to the species richness of the community with the highest species richness. Set PSVmncd=NULL, PSVpool=NULL and run pcd again.")
  }

	# m=number of communities; n=number of species; nsr=maximum sr value across all communities
	m <- dim(comm)[1]
	n <- dim(comm)[2]
  nsr <-max(rowSums(comm))
  if(is.null(SSii) & is.null(SCii))   #If the user already has calculated the mean conditional PSV values for all levels of SR
  {                                   #and the PSV of the species pool
  	SSii <- array(0,nsr)
  	n1 <- 2
  	for (n2 in 1:nsr)
  	{
  		temp <- array(0,reps)
  		for (t in 1:reps)
  		{
  			rp <- sample(n)
  			pick1 <- rp[1:n1]

  			rp <- sample(n)
  			pick2 <- rp[1:n2]

  			C11 <- V[pick1,pick1]
  			C22 <- V[pick2,pick2]
  			C12 <- V[pick1,pick2]

  			invC22 <- solve(C22)
  			S11 <- C11 - C12%*%invC22%*%t(C12)
  			SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1-1))
  			temp[t] <- SS11
  		}
  		SSii[n2] <- mean(temp)
  	}
  	SCii=1-(sum(V)-sum(diag(V)))/(n*(n-1))
  }

	# calculate PCD
	PCD <- array(NA,c(m,m))
	PCDc <- array(NA,c(m,m))
	PCDp <- array(NA,c(m,m))
	for (i in 1:(m-1))
	{
		for (j in (i+1):m)
		{
			pick1 <- (1:n)[comm[i,]==1]
 			pick2 <- (1:n)[comm[j,]==1]

			n1 <- length(pick1)
			n2 <- length(pick2)

			C <- V[c(pick1, pick2),c(pick1, pick2)]

			C11 <- C[1:n1,1:n1]
			C22 <- C[(n1+1):(n1+n2),(n1+1):(n1+n2)]
			C12 <- C[1:n1,(n1+1):(n1+n2)]
			if(is.null(dim(C12)))
      {
        if(is.null(dim(C22))){C12<-as.matrix(C12)} else {C12<-t(as.matrix(C12))}
      }

			invC11 <- solve(C11)
			S22 <- C22 - t(C12)%*%invC11%*%C12

			invC22 <- solve(C22)
			S11 <- C11 - C12%*%invC22%*%t(C12)
      if(n1>1)
      {
			 SC11 <- (n1*sum(diag(C11))-sum(C11))/(n1*(n1-1))
       SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1-1))
			} else {
       SC11 <- (n1*sum(diag(C11))-sum(C11))/(n1*(n1))
       SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1))
      }
      if(n2>1)
      {
        SC22 <- (n2*sum(diag(C22))-sum(C22))/(n2*(n2-1))
        SS22 <- (n2*sum(diag(S22))-sum(S22))/(n2*(n2-1))
      } else {
        SC22 <- (n2*sum(diag(C22))-sum(C22))/(n2*(n2))
        SS22 <- (n2*sum(diag(S22))-sum(S22))/(n2*(n2))
      }

			D=(n1*SS11 + n2*SS22)/(n1*SC11 + n2*SC22)

			a <- length(unique(c(pick1, pick2)))
			b <- length(pick1)-a
			cc <- length(pick2)-a
			dsor <- 2*a/(2*a+b+cc) - 1

			pred.D <- (n1*SSii[n2]+n2*SSii[n1])/(n1*SCii+n2*SCii)
			pred.dsor <- 1 - 2*n1*n2/((n1+n2)*n)

			PCD[i,j] <- D/pred.D
			PCDc[i,j] <- dsor/pred.dsor
			PCDp[i,j] <- PCD[i,j]/PCDc[i,j]
		}
	}
	colnames(PCD)<-rownames(comm)
  rownames(PCD)<-rownames(comm)
  colnames(PCDc)<-rownames(comm)
  rownames(PCDc)<-rownames(comm)
  colnames(PCDp)<-rownames(comm)
  rownames(PCDp)<-rownames(comm)

  return(list(PCD=as.dist(t(PCD)), PCDc=as.dist(t(PCDc)), PCDp=as.dist(t(PCDp)),  PSVmncd=SSii, PSVpool=SCii))
}
