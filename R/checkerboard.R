`checkerboard` <-
function(x) {
	#Gotelli 2000: Checker = Sum (Si - Q)(Sk - Q) / ((R*(R-1))/2)
	#where Si = total for row(species) i, R = num rows(spp), Q = num sites where both spp present
	x <- decostand(x,method="pa")
	Nsites <- dim(x)[1]
	S <- apply(x,2,sum)
	R <- length(S)
	Checker.ij <- matrix(nrow=R,ncol=R,dimnames=list(colnames(x),colnames(x)))
	for (i in 1:R) {
		for (j in 1:R) {
			Q <- sum(x[,i]*x[,j])
			Checker.ij[i,j] <- ((S[i] - Q)*(S[j] - Q)) / ((R*(R-1))/2)
		}
	}
	as.dist(Checker.ij)
}

