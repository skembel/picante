raoD <- function (comm, phy=NULL)
{
    res <- list()

	if (is.null(phy)) {
	    tij <- 1-diag(x=rep(1,length(comm[1,])))
	} else {
	    if (!is.ultrametric(phy)) {
	        stop("Phylogeny must be ultrametric")
	    }
	    dat <- match.phylo.comm(phy, comm)
	    comm <- dat$comm
	    phy <- dat$phy
	    tij <- cophenetic(phy) / 2
	}

    x <- as.matrix(comm)    
    S <- length(x[1,])
    N <- length(x[,1])
    total <- apply(x, 1, sum)
    samp.relabund <- total / sum(x)
    x.combined <- matrix(apply(x, 2, sum), nrow=1) / sum(x)
    x <- sweep(x, 1, total, "/")

	D <- vector(length=N)
    names(D) <- rownames(x)
	for (k in 1:N)
		D[k] <-
			sum ( tij * outer(as.vector(t(x[k,])),as.vector(t(x[k,]))) )
    res$Dkk <- D

	Dkl <- matrix(nrow=N,ncol=N)    
	for (k in 1:N) {
		for (l in 1:N) {
			Dkl[k,l] <-
			sum ( tij * outer(as.vector(t(x[k,])),as.vector(t(x[l,]))) )
		}
	}

    row.names(Dkl) <- row.names(x)
    colnames(Dkl) <- row.names(x)
	H <- Dkl
	res$Dkl <- Dkl

	for (k in 1:N) {
		for (l in 1:N) {
			H[k,l] <- Dkl[k,l] - (Dkl[k,k] + Dkl[l,l]) / 2
		}
	}
    res$H <- H
    
    res$total <- sum ( tij * outer(as.vector(t(x.combined)),as.vector(t(x.combined))) )

    res$alpha <- sum(res$Dkk * samp.relabund)

    res$beta <- res$total - res$alpha
    
    res$Fst <- res$beta / res$total
    
    return(res)
    
}
