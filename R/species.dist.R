`species.dist` <-
function (x, metric=c("cij","jaccard","checkerboard","doij")) {
	metric <- match.arg(metric)
    if (identical(metric,"checkerboard")) {
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
        return(as.dist(Checker.ij))
    }
    if (identical(metric,"cij")) {
        #Schoener index of co-occurrence
        x <- decostand(x,method="total",MARGIN=2)
        return(1 - (0.5 * dist(t(x),method="manhattan")))
    }    
    if (identical(metric,"jaccard")) {    
        return( 1 - vegdist(t(sortColumns(x)), method = "jaccard"))
    }
    if (identical(metric,"doij")) {
        #Hardy's standardized version of checkerboard
        #doij = (Pij - Pi*Pj)/(Pi*Pj)
        x <- as.matrix(decostand(x,method="pa"))
        Nsites <- dim(x)[1]
        P <- apply(x,2,sum) / Nsites
        N <- length(P)
        doij <- matrix(nrow=N,ncol=N,dimnames=list(colnames(x),colnames(x)))
        for (i in 1:N-1) {
            for (j in (i+1):N) {
                Pij <- sum(x[,i]*x[,j])/Nsites
                doij[i,j] <- ((Pij - (P[i]*P[j]))/(P[i]*P[j]))
            }
        }
        return(as.dist(t(doij)))
    }
}
