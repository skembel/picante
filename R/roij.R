`roij` <-
function(x) {
	#Hardy's standardized version of checkerboard
	#roij = (Pij - Pi*Pj)/(Pi*Pj)
	x <- as.matrix(decostand(x,method="pa"))
	Nsites <- dim(x)[1]
	P <- apply(x,2,sum) / Nsites
	N <- length(P)
	roij <- matrix(nrow=N,ncol=N,dimnames=list(colnames(x),colnames(x)))
	for (i in 1:N-1) {
		for (j in (i+1):N) {
			Pij <- sum(x[,i]*x[,j])/Nsites
			roij[i,j] <- ((Pij - (P[i]*P[j]))/(P[i]*P[j]))
		}
	}
	as.dist(t(roij))
}

